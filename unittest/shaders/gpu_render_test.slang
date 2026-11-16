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
    uint*   cmds;           // draw commands a kernel writes (gen_cmd)
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

// Conservative depth, as md_render writes it (reversed Z: nearer is larger):
// SV_DepthLessEqual, a depth no larger than the rasterised one, here that
// scaled by color.x (<= 1). With flags & 8, odd columns are discarded.
struct DepthLessOut {
    float4 color : SV_Target0;
    float  depth : SV_DepthLessEqual;
};

[shader("fragment")]
DepthLessOut depth_less_fs(VOut i) {
    RArgs a = MD_ARGS;
    if ((a.flags & 8u) != 0u && (uint(i.pos.x) & 1u) != 0u) discard;
    DepthLessOut o;
    o.color = i.col;
    o.depth = i.pos.z * a.color.x;
    return o;
}

// Early depth rejection (gpu_render.early_depth_rejection): full-screen
// layers whose fragments run a long loop (a.width turns), at depth
// a.color.y + a.color.z * instance. They write a.id + instance, so the
// image shows which layer won, and the loop's result decides nothing but
// cannot be dropped. One fragment stage per way of ending: as md_render's
// spline and polygon fragments end (a depth no nearer than the box's, and a
// discard where the ray misses) and the ways around it.
struct LOut {
    float4 pos : SV_Position;
    nointerpolation uint inst : INST;
};

[shader("vertex")]
LOut layer_vs(uint vid : SV_VertexID, uint iid : SV_InstanceID) {
    RArgs a = MD_ARGS;
    LOut  o;
    o.pos  = float4(vid == 1u ? 3.0 : -1.0, vid == 2u ? 3.0 : -1.0, a.color.y + a.color.z * float(iid), 1.0);
    o.inst = iid;
    return o;
}

uint layer_work(LOut i) {
    RArgs a = MD_ARGS;
    float v = i.pos.x * 0.37 + i.pos.y * 0.11;
    for (uint k = 0u; k < a.width; ++k) v = sin(v) * 1.3 + cos(v * 0.7) + 0.01;
    return a.id + i.inst + (asuint(v) == 0x7FC00001u ? 1u : 0u);   // (never: a NaN of its own)
}

[shader("fragment")]
uint layer_plain_fs(LOut i) : SV_Target { return layer_work(i); }

[shader("fragment")]
uint layer_discard_fs(LOut i) : SV_Target {
    uint id = layer_work(i);
    if (id == 0xFFFFFFFFu) discard;
    return id;
}

struct LDepthAny  { uint id : SV_Target; float depth : SV_Depth; };
struct LDepthLess { uint id : SV_Target; float depth : SV_DepthLessEqual; };
struct LDepthMore { uint id : SV_Target; float depth : SV_DepthGreaterEqual; };

[shader("fragment")]
LDepthAny layer_depth_any_fs(LOut i) {
    LDepthAny o; o.id = layer_work(i); o.depth = i.pos.z; return o;
}

// reversed Z (GREATER_EQUAL): no nearer is smaller
[shader("fragment")]
LDepthLess layer_depth_less_fs(LOut i) {
    LDepthLess o; o.id = layer_work(i); o.depth = i.pos.z * 0.999; return o;
}

[shader("fragment")]
LDepthLess layer_depth_less_discard_fs(LOut i) {
    LDepthLess o; o.id = layer_work(i);
    if (o.id == 0xFFFFFFFFu) discard;
    o.depth = i.pos.z * 0.999;
    return o;
}

// Forced early tests: the depth test before the stage, which then may not
// write depth (the pipeline does not), and a side effect (a.cmds: an
// atomic max a pixel, 1024 slots), as a visibility kept by atomics would.
[earlydepthstencil]
[shader("fragment")]
uint layer_early_fs(LOut i) : SV_Target {
    RArgs a = MD_ARGS;
    uint  id = layer_work(i);
    InterlockedMax(a.cmds[(uint(i.pos.y) * 512u + uint(i.pos.x)) & 1023u], id);
    return id;
}

// standard Z (LESS_EQUAL): no nearer is larger
[shader("fragment")]
LDepthMore layer_depth_greater_fs(LOut i) {
    LDepthMore o; o.id = layer_work(i); o.depth = min(i.pos.z * 1.001 + 1e-6, 1.0); return o;
}

// Rendered texture -> buffer, through the bindless heap.
[shader("compute")]
[numthreads(8, 8, 1)]
void read_tex(uint3 tid : SV_DispatchThreadID) {
    RArgs a = MD_ARGS;
    if (tid.x >= a.width || tid.y >= a.height) return;
    a.colors[tid.y * a.width + tid.x] = a.tex.Load(int3(int2(tid.xy), 0));
}

// Draw commands produced by a kernel and consumed by an indirect draw in the
// same stream, as md_render's culls write its rasters': a draw of 3 vertices
// at cmds[0], an indexed draw of 3 indices at cmds[4] (md_gpu_draw_cmd_t,
// md_gpu_draw_indexed_cmd_t), one instance each.
[shader("compute")]
[numthreads(1, 1, 1)]
void gen_cmd(uint3 tid : SV_DispatchThreadID) {
    RArgs a = MD_ARGS;
    uint  n = (a.flags & 16u) != 0u ? min(a.cmds[32], 1u) : 1u;   // 16: the instance count gen_cnt counted
    a.cmds[0] = 3u; a.cmds[1] = n; a.cmds[2] = 0u; a.cmds[3] = 0u;
    a.cmds[4] = 3u; a.cmds[5] = n; a.cmds[6] = 0u; a.cmds[7] = 0u; a.cmds[8] = 0u;
}

// A count made with atomics, which gen_cmd turns into a command (md_render's
// cull and finalize)
[shader("compute")]
[numthreads(64, 1, 1)]
void gen_cnt(uint3 tid : SV_DispatchThreadID) {
    RArgs a = MD_ARGS;
    if (tid.x == 0u) InterlockedAdd(a.cmds[32], 1u);
}

// Positions produced by a kernel and consumed by a draw in the same stream.
[shader("compute")]
[numthreads(64, 1, 1)]
void gen_pos(uint3 tid : SV_DispatchThreadID) {
    RArgs a = MD_ARGS;
    if (tid.x >= a.width) return;
    a.pos[tid.x] = a.colors[tid.x] * a.color;
}
