// gpu_render_test.slang -- vertex, fragment and compute entries used by the
// rendering tests in unittest/test_gpu.c.
//
// Every entry shares one argument struct, as a draw's two stages must. The
// vertex shaders pull clip-space positions through a pointer; there are no
// vertex buffers anywhere in md_gpu.

#include "../../src/shaders/md_gpu.slang"

struct Draw {
    float4 offset;          // xy added to the position
    float4 color;           // multiplies the colour
};

struct RArgs {
    float4  color;          // constant colour; .x is the depth written by depth_fs
    float4* pos;            // per-vertex clip-space positions
    float4* colors;         // per-vertex colours (flags & 1)
    Draw*   draws;          // per-instance records (flags & 2), per-draw when indirect
    uint    id;             // base of the value written to integer targets
    uint    flags;          // 1: per-vertex colour, 2: per-instance record, 4: sample tex
    DescriptorHandle<Texture2D>    tex;
    DescriptorHandle<SamplerState> smp;
    uint    width;
    uint    height;
    uint    _p0;
    uint    _p1;
};
MD_SHADER_ARGS(RArgs);

struct VOut {
    float4 pos : SV_Position;
    float4 col : COLOR;
    float2 uv  : UV;
    nointerpolation uint inst : INST;
};

VOut shade(uint v, uint rec) {
    RArgs  a = MD_ARGS;
    VOut   o;
    float4 p = a.pos[v];
    float4 c = a.color;
    if ((a.flags & 1u) != 0u) c = a.colors[v];
    if ((a.flags & 2u) != 0u) {
        Draw d = a.draws[rec];
        p.xy += d.offset.xy;
        c *= d.color;
    }
    o.pos  = p;
    o.col  = c;
    o.uv   = float2(p.x * 0.5 + 0.5, 0.5 - p.y * 0.5);   // top-left origin, like the framebuffer
    o.inst = rec;
    return o;
}

// Direct draws: SV_VertexID / SV_InstanceID are draw-relative everywhere.
[shader("vertex")]
VOut basic_vs(uint vid : SV_VertexID, uint iid : SV_InstanceID) {
    return shade(vid, iid);
}

// Indirect draws: first_vertex picks the vertices, first_instance the record.
[shader("vertex")]
VOut indirect_vs(uint vid : SV_VertexID, uint iid : SV_InstanceID,
                 uint base_vertex : SV_StartVertexLocation, uint base_instance : SV_StartInstanceLocation) {
    uint v = base_vertex + md_draw_vertex(vid, base_vertex);
    uint i = md_draw_instance(iid, base_instance);
    VOut o = shade(v, base_instance);
    o.inst = base_instance * 16u + i;
    return o;
}

struct PVOut {
    float4 pos : SV_Position;
    float4 col : COLOR;
    float2 uv  : UV;
    nointerpolation uint inst : INST;
    float  size : SV_PointSize;
};

[shader("vertex")]
PVOut point_vs(uint vid : SV_VertexID) {
    VOut  v = shade(vid, 0u);
    PVOut o;
    o.pos  = v.pos;
    o.col  = v.col;
    o.uv   = v.uv;
    o.inst = v.inst;
    o.size = 1.0;
    return o;
}

[shader("fragment")]
float4 color_fs(VOut i) : SV_Target0 {
    RArgs a = MD_ARGS;
    float4 c = i.col;
    if ((a.flags & 4u) != 0u) c *= a.tex.SampleLevel(a.smp, i.uv, 0.0);
    return c;
}

struct MrtOut {
    float4 color : SV_Target0;
    uint   id    : SV_Target1;
};

[shader("fragment")]
MrtOut mrt_fs(VOut i) {
    RArgs  a = MD_ARGS;
    MrtOut o;
    o.color = i.col;
    o.id    = a.id + i.inst;
    return o;
}

[shader("fragment")]
uint id_fs(VOut i) : SV_Target0 {
    return MD_ARGS.id + i.inst;
}

struct DepthOut {
    float4 color : SV_Target0;
    float  depth : SV_Depth;
};

[shader("fragment")]
DepthOut depth_fs(VOut i) {
    DepthOut o;
    o.color = i.col;
    o.depth = MD_ARGS.color.x;
    return o;
}

// Rendered texture -> buffer, through the bindless heap.
[shader("compute")]
[numthreads(8, 8, 1)]
void read_tex(uint3 tid : SV_DispatchThreadID) {
    RArgs a = MD_ARGS;
    if (tid.x >= a.width || tid.y >= a.height) return;
    a.colors[tid.y * a.width + tid.x] = a.tex.Load(int3(int2(tid.xy), 0));
}

// Positions produced by a kernel and consumed by a draw in the same stream.
[shader("compute")]
[numthreads(64, 1, 1)]
void gen_pos(uint3 tid : SV_DispatchThreadID) {
    RArgs a = MD_ARGS;
    if (tid.x >= a.width) return;
    a.pos[tid.x] = a.colors[tid.x] * a.color;
}
