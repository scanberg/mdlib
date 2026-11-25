// md_post_overlay.slang -- md_post_overlay's (md_post.h): the caller's marks over
// the image shown, a quad (two triangles) per item, alpha-blended in order.
// Discs and segments are implicit: a fragment's distance to the disc's centre
// or to the segment, in pixels, gives its coverage and its colour (the fill,
// the outline), so their edges are smooth at any size. Triangles are drawn as
// rasterised.
//
// Pixels are x right, y down, their centres at + 0.5 (md_gpu's framebuffer).
// Out premultiplied (blended ONE, ONE_MINUS_SRC_ALPHA). Layout rules of
// md_gpu.slang. Mirrored in md_post.c.

#include "../md_gpu.slang"

struct Item {                   // md_post_overlay_item_t
    float4 p0;                  // xyz; w: the kind (its bits)
    float4 p1;                  // xyz; w: the flags (its bits)
    float4 p2;                  // xyz; w: the size
    uint   color, color2;
    float  width;
    uint   _pad;
};

struct Args {
    float4 v0, v1, v2, v3;      // world to view (columns)
    float4 q0, q1, q2, q3;      // view to clip (columns)
    float4 sz;                  // the target's width, height; the depth's width, height
    float4 nf;                  // x: the near plane's distance in view (perspective), y: 1 orthographic
    Item*  items;
    uint   count, has_depth;
    DescriptorHandle<Texture2D<float>> depth;
    uint   _p0, _p1;
};
MD_SHADER_ARGS(Args);

static const uint KIND_DISC     = 0u;
static const uint KIND_SEGMENT  = 1u;
static const uint KIND_TRIANGLE = 2u;
static const uint F_PIXELS = 1u << 0;
static const uint F_DEPTH  = 1u << 1;

float4 unpack_rgba8(uint c) {
    return float4(float(c & 0xFFu), float((c >> 8) & 0xFFu), float((c >> 16) & 0xFFu), float(c >> 24)) / 255.0;
}

float3 to_view(Args a, float3 p) { return (a.v0 * p.x + a.v1 * p.y + a.v2 * p.z + a.v3).xyz; }
float4 to_clip(Args a, float3 v) { return a.q0 * v.x + a.q1 * v.y + a.q2 * v.z + a.q3; }

// A clip position's pixel (x right, y down)
float2 to_px(Args a, float4 c) {
    float2 n = c.xy / c.w;
    return float2((n.x * 0.5 + 0.5) * a.sz.x, (0.5 - n.y * 0.5) * a.sz.y);
}
float4 px_clip(Args a, float2 px) {
    return float4(px.x / a.sz.x * 2.0 - 1.0, 1.0 - px.y / a.sz.y * 2.0, 0.5, 1.0);
}

// Pixels per world unit at a clip w (view depth, or 1 orthographic)
float px_per_world(Args a, float w) { return a.q1.y * 0.5 * a.sz.y / w; }

struct OV {
    float4 pos : SV_Position;
    nointerpolation float4 g0 : TEXCOORD0;   // disc: its centre (px), radius, outline width (px); segment: its ends (px)
    nointerpolation float4 g1 : TEXCOORD1;   // segment: its half width at each end (px), its ends' view z
    nointerpolation uint   it : TEXCOORD2;   // the item; past the count: culled
    float  z : TEXCOORD3;                    // a triangle's view z (a disc's: g1.z)
};

static const uint QUAD[6] = { 0u, 1u, 2u, 0u, 2u, 3u };

OV culled() {
    OV o;
    o.pos = float4(0.0, 0.0, -2.0, 1.0);
    o.g0 = o.g1 = float4(0.0);
    o.it = 0xFFFFFFFFu;
    o.z = 0.0;
    return o;
}

// The near plane's: a point of view z in front of it moved along to it
float3 clip_near(Args a, float3 p, float3 q) {
    float zn = -a.nf.x * 1.001;
    if (a.nf.y > 0.5 || p.z <= zn) return p;
    return p + (q - p) * ((zn - p.z) / (q.z - p.z));
}

[shader("vertex")]
OV overlay_vs(uint vid : SV_VertexID) {
    Args a = MD_ARGS;
    uint i = vid / 6u, k = vid - i * 6u;
    if (i >= a.count) return culled();
    Item it = a.items[i];
    uint kind = asuint(it.p0.w), flags = asuint(it.p1.w);
    float size = it.p2.w;
    OV o = culled();
    o.it = i;

    if (kind == KIND_TRIANGLE) {
        float3 v = to_view(a, k == 1u ? it.p1.xyz : (k == 2u ? it.p2.xyz : it.p0.xyz));
        o.pos = k < 3u ? to_clip(a, v) : float4(0.0, 0.0, -2.0, 1.0);
        o.z = v.z;
        return o;
    }

    float2 corner = float2(QUAD[k] == 1u || QUAD[k] == 2u ? 1.0 : -1.0, QUAD[k] >= 2u ? 1.0 : -1.0);
    if (kind == KIND_DISC) {
        float3 c = to_view(a, it.p0.xyz);
        if (a.nf.y < 0.5 && c.z > -a.nf.x) return culled();
        float4 cc = to_clip(a, c);
        float  s  = (flags & F_PIXELS) != 0u ? 1.0 : px_per_world(a, cc.w);
        float  r  = size * s, w = it.width * s;
        float2 cp = to_px(a, cc);
        o.pos = px_clip(a, cp + corner * (r + 1.0));
        o.g0 = float4(cp, r, w);
        o.g1 = float4(0.0, 0.0, c.z, 0.0);
        return o;
    }

    // a segment: its ends in front of the near plane, a box about them in pixels
    float3 va = to_view(a, it.p0.xyz), vb = to_view(a, it.p1.xyz);
    if (a.nf.y < 0.5 && va.z > -a.nf.x && vb.z > -a.nf.x) return culled();
    float3 ca = clip_near(a, va, vb), cb = clip_near(a, vb, va);
    float4 qa = to_clip(a, ca), qb = to_clip(a, cb);
    float2 pa = to_px(a, qa), pb = to_px(a, qb);
    float  ra = (flags & F_PIXELS) != 0u ? size : size * px_per_world(a, qa.w);
    float  rb = (flags & F_PIXELS) != 0u ? size : size * px_per_world(a, qb.w);
    float2 d  = pb - pa;
    float  len = length(d);
    float2 u  = len > 1e-6 ? d / len : float2(1.0, 0.0);
    float2 n  = float2(-u.y, u.x);
    float  e  = max(ra, rb) + 1.0;
    float2 p  = (corner.x < 0.0 ? pa - u * e : pb + u * e) + n * (corner.y * e);
    o.pos = px_clip(a, p);
    o.g0 = float4(pa, pb);
    o.g1 = float4(ra, rb, ca.z, cb.z);
    return o;
}

// Whether the scene is in front of view z at the pixel
bool behind(Args a, float2 px, float z) {
    int2  q = int2(px * a.sz.zw / a.sz.xy);
    q = clamp(q, int2(0), int2(a.sz.zw) - 1);
    float d = a.depth.Load(int3(q, 0));
    if (d <= 0.0) return false;
    // the view z of reverse-Z depth d: d (P23 z + P33) = P22 z + P32
    float zs = (a.q3.z - d * a.q3.w) / (d * a.q2.w - a.q2.z);
    return zs > z + 1e-3 * abs(z) + 1e-4;
}

[shader("fragment")]
float4 overlay_fs(OV i) : SV_Target0 {
    Args a = MD_ARGS;
    if (i.it >= a.count) discard;
    Item it = a.items[i.it];
    uint kind = asuint(it.p0.w), flags = asuint(it.p1.w);
    float2 px = i.pos.xy;
    float4 c = unpack_rgba8(it.color);
    float  z;
    float4 o;
    if (kind == KIND_TRIANGLE) {
        z = i.z;
        o = float4(c.rgb * c.a, c.a);
    } else if (kind == KIND_DISC) {
        z = i.g1.z;
        float r = i.g0.z, w = i.g0.w;
        float d = length(px - i.g0.xy);
        float cov  = saturate(r - d + 0.5);
        float fill = w > 0.0 ? saturate(r - w - d + 0.5) : cov;
        float4 c2 = unpack_rgba8(it.color2);
        float  af = c.a * fill;
        float  ao = w > 0.0 ? c2.a * cov * (1.0 - af) : 0.0;
        o = float4(c.rgb * af + c2.rgb * ao, af + ao);
    } else {
        float2 pa = i.g0.xy, pb = i.g0.zw, d = pb - pa;
        float  l2 = dot(d, d);
        float  t  = l2 > 1e-12 ? saturate(dot(px - pa, d) / l2) : 0.0;
        float  r  = lerp(i.g1.x, i.g1.y, t);
        float  cov = saturate(r - length(px - (pa + d * t)) + 0.5);
        z = lerp(i.g1.z, i.g1.w, t);
        o = float4(c.rgb, 1.0) * (c.a * cov);
    }
    if (o.a <= 0.0) discard;
    if ((flags & F_DEPTH) != 0u && a.has_depth != 0u && behind(a, px, z)) discard;
    return o;
}
