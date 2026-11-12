// md_render.slang -- raster stages and kernels of md_render.c.
//
// Raster writes depth and a 64-bit id per pixel; the resolve, mask and pick
// kernels intersect the winning primitive again with the very same functions
// the fragment stages used. That shared code is what a visibility buffer rests
// on: raster and resolve must agree on which primitive and where.
//
// Every struct here is mirrored in md_render.c. Layout rules of md_gpu.slang
// apply: no 3-vectors, 4-vectors at 16-byte offsets, pointers at 8.

#include "../md_gpu.slang"

// ---------------------------------------------------------------------------
// Shared data
// ---------------------------------------------------------------------------

struct Mat {                    // column-major 4x4
    float4 c0;
    float4 c1;
    float4 c2;
    float4 c3;
};

float4 mul4(Mat m, float4 v) { return m.c0 * v.x + m.c1 * v.y + m.c2 * v.z + m.c3 * v.w; }
float3 mulp(Mat m, float3 p) { return (m.c0.xyz * p.x + m.c1.xyz * p.y + m.c2.xyz * p.z + m.c3.xyz); }
float3 mulv(Mat m, float3 v) { return (m.c0.xyz * v.x + m.c1.xyz * v.y + m.c2.xyz * v.z); }

// Matrices as their four columns: the argument-layout check models vectors,
// not nested structs. M(v, name) rebuilds one.
struct ViewData {
    float4 view0, view1, view2, view3;                  // world -> view
    float4 proj0, proj1, proj2, proj3;                  // view -> clip, as rasterised (jittered)
    float4 inv_proj0, inv_proj1, inv_proj2, inv_proj3;  // of proj
    float4 inv_view0, inv_view1, inv_view2, inv_view3;  // view -> world
    float4 view_prev0, view_prev1, view_prev2, view_prev3;// of the previous raster
    float4 proj_prev0, proj_prev1, proj_prev2, proj_prev3;// of the previous raster, unjittered
    float4 proj_cur0, proj_cur1, proj_cur2, proj_cur3;  // unjittered
    float4 viewport;            // w, h, 1/w, 1/h
    float  dt;                  // time since the previous raster
    uint   ortho;
    uint   batch_count;
    uint   _pad;
};

Mat mat_of(float4 a, float4 b, float4 c, float4 d) { Mat m; m.c0 = a; m.c1 = b; m.c2 = c; m.c3 = d; return m; }
#define M(v, n) mat_of(v.n##0, v.n##1, v.n##2, v.n##3)

// Batch.flags
static const uint F_INDEX        = 1u << 0;
static const uint F_VELOCITY     = 1u << 1;
static const uint F_RADIUS       = 1u << 2;
static const uint F_ATTR_INDEX   = 1u << 3;
static const uint F_VISIBLE0     = 1u << 4;
static const uint F_VISIBLE1     = 1u << 5;
static const uint F_STATE0       = 1u << 6;
static const uint F_STATE1       = 1u << 7;
static const uint F_ID           = 1u << 8;
static const uint F_RADIUS_SCALE = 1u << 9;
static const uint F_COLOR        = 1u << 10;
static const uint F_PALETTE_IDX  = 1u << 11;

// md_render_kind_t
static const uint KIND_SPHERES  = 1u;
static const uint KIND_CAPSULES = 2u;
static const uint KIND_EXTERNAL = 6u;

// md_render_color_mode_t
static const uint COLOR_POINTS  = 0u;
static const uint COLOR_UNIFORM = 1u;
static const uint COLOR_PALETTE = 2u;

struct Batch {
    uint    kind;
    uint    count;
    uint    instance_count;
    uint    flags;

    uint    stride;             // in bytes
    uint    vel_stride;
    uint    attr_count;
    uint    color_mode;

    uint*   position;
    uint*   velocity;
    float*  radius;
    uint*   attr_index;
    uint*   color;
    uint*   palette_index;      // uint16 per attribute, two per word
    uint*   palette;
    uint*   visible0;
    uint*   visible1;
    uint*   state0;
    uint*   state1;
    uint*   index;
    uint*   id;
    float*  radius_scale;
    float4* transforms;         // 4 columns per instance
    float4* transforms_prev;

    float   s0;                 // SPHERES: radius_scale  CAPSULES: radius
    float   s1;                 // SPHERES: radius_add
    uint    uniform_color;
    uint    palette_count;

    float4  tint;

    float   sharpness;
    float   desaturate;
    uint    _pad0;
    uint    _pad1;
};

// One argument struct for every entry point (md_gpu shares one root between a
// draw's two stages, and keeping kernels on it too means one C mirror).
//   raster:  view, batches, slot (the batch's position in the frame), mode
//            (0 visibility, 1 mask)
//   kernels: view, batches, the targets, size; pick: px, py, pick_out
struct Args {
    ViewData* view;
    Batch*    batches;
    uint      slot;
    uint      mode;
    DescriptorHandle<Texture2D<uint2>>     id;
    DescriptorHandle<RWTexture2D<float4>>  color;
    DescriptorHandle<RWTexture2D<float2>>  normal;
    DescriptorHandle<RWTexture2D<float2>>  motion;
    DescriptorHandle<RWTexture2D<float4>>  mask;
    uint      width;
    uint      height;
    uint      px;
    uint      py;
    uint*     pick_out;
};
MD_SHADER_ARGS(Args);

// ---------------------------------------------------------------------------
// Fetch
// ---------------------------------------------------------------------------

float3 load_pos(Batch b, uint i) {
    uint* p = b.position + i * (b.stride >> 2);
    return float3(asfloat(p[0]), asfloat(p[1]), asfloat(p[2]));
}

float3 load_vel(Batch b, uint i) {
    if ((b.flags & F_VELOCITY) == 0u) return float3(0.0);
    uint* p = b.velocity + i * (b.vel_stride >> 2);
    uint w0 = p[0];
    uint w1 = p[1];
    return float3(f16tof32(w0 & 0xFFFFu), f16tof32(w0 >> 16), f16tof32(w1 & 0xFFFFu));
}

uint attr_of(Batch b, uint pi) {
    return ((b.flags & F_ATTR_INDEX) != 0u) ? b.attr_index[pi] : pi;
}

bool test_bit(uint* set, uint i) {
    return ((set[i >> 5] >> (i & 31u)) & 1u) != 0u;
}

bool attr_visible(Batch b, uint a) {
    if ((b.flags & F_VISIBLE0) != 0u && !test_bit(b.visible0, a)) return false;
    if ((b.flags & F_VISIBLE1) != 0u && !test_bit(b.visible1, a)) return false;
    return true;
}

uint attr_state(Batch b, uint a) {
    uint s = 0u;
    if ((b.flags & F_STATE0) != 0u && test_bit(b.state0, a)) s |= 1u;
    if ((b.flags & F_STATE1) != 0u && test_bit(b.state1, a)) s |= 2u;
    return s;
}

// Point k (0 or 1) of element e with `per` points per element
uint point_of(Batch b, uint e, uint k, uint per) {
    uint i = e * per + k;
    return ((b.flags & F_INDEX) != 0u) ? b.index[i] : i;
}

Mat load_transform(float4* t, uint inst) {
    Mat m;
    m.c0 = t[inst * 4u + 0u];
    m.c1 = t[inst * 4u + 1u];
    m.c2 = t[inst * 4u + 2u];
    m.c3 = t[inst * 4u + 3u];
    return m;
}

float sphere_radius(Batch b, uint pi) {
    float r = ((b.flags & F_RADIUS) != 0u) ? b.radius[pi] : 0.0;
    return r * b.s0 + b.s1;
}

float capsule_radius(Batch b, uint e) {
    float r = b.s0;
    if ((b.flags & F_RADIUS_SCALE) != 0u) r *= b.radius_scale[e];
    return r;
}

// ---------------------------------------------------------------------------
// Rays and intersections (view space)
// ---------------------------------------------------------------------------

struct Ray { float3 o; float3 d; };

float3 unproject(Mat inv_proj, float2 ndc, float z) {
    float4 p = mul4(inv_proj, float4(ndc, z, 1.0));
    return p.xyz / p.w;
}

// The ray of pixel position `pix` (pixel units, top-left origin; a pixel's
// centre is at +0.5), through the projection as it was rasterised.
Ray pixel_ray(ViewData v, float2 pix) {
    float2 ndc = float2(pix.x * v.viewport.z * 2.0 - 1.0, 1.0 - pix.y * v.viewport.w * 2.0);
    Ray r;
    r.o = unproject(M(v, inv_proj), ndc, 1.0);          // near plane under reverse Z
    float3 f = unproject(M(v, inv_proj), ndc, 0.5);
    r.d = normalize(f - r.o);
    return r;
}

// Nearest t >= 0, or -1. `h_clamped`: return the closest approach for a miss
// (resolve must not lose a pixel the raster kept to a rounding difference).
float sphere_hit(Ray r, float3 c, float rad, bool h_clamped) {
    float3 oc = r.o - c;
    float  b  = dot(oc, r.d);
    float  cc = dot(oc, oc) - rad * rad;
    float  h  = b * b - cc;
    if (h < 0.0) {
        if (!h_clamped) return -1.0;
        h = 0.0;
    }
    float t = -b - sqrt(h);
    return t;
}

// Ray against the capsule pa-pb of radius rad (Quilez). t, or -1.
float capsule_hit(Ray r, float3 pa, float3 pb, float rad) {
    float3 ba   = pb - pa;
    float3 oa   = r.o - pa;
    float  baba = dot(ba, ba);
    float  bard = dot(ba, r.d);
    float  baoa = dot(ba, oa);
    float  rdoa = dot(r.d, oa);
    float  oaoa = dot(oa, oa);
    float  a    = baba - bard * bard;
    float  b    = baba * rdoa - baoa * bard;
    float  c    = baba * oaoa - baoa * baoa - rad * rad * baba;
    float  h    = b * b - a * c;
    if (h >= 0.0 && a > 1e-12 * baba) {
        float t = (-b - sqrt(h)) / a;
        float y = baoa + t * bard;
        if (y > 0.0 && y < baba) return t;
    }
    // the caps
    float best = -1.0;
    float t0 = sphere_hit(r, pa, rad, false);
    if (t0 >= 0.0) best = t0;
    float t1 = sphere_hit(r, pb, rad, false);
    if (t1 >= 0.0 && (best < 0.0 || t1 < best)) best = t1;
    return best;
}

float capsule_param(float3 p, float3 pa, float3 pb) {
    float3 ba = pb - pa;
    return saturate(dot(p - pa, ba) / max(dot(ba, ba), 1e-20));
}

float clip_depth(ViewData v, float3 p) {
    float4 c = mul4(M(v, proj), float4(p, 1.0));
    return c.z / c.w;
}

// ---------------------------------------------------------------------------
// Ids
// ---------------------------------------------------------------------------

uint2 encode_id(uint slot, uint part, uint64_t linear) {
    uint hi = uint(linear >> 32) & 0x3FFFFu;
    return uint2(uint(linear & 0xFFFFFFFFull), ((slot + 1u) << 20) | ((part & 3u) << 18) | hi);
}

// ---------------------------------------------------------------------------
// Raster
// ---------------------------------------------------------------------------

static const float4 CULLED = float4(0.0, 0.0, -2.0, 1.0);   // outside the depth range: clipped

struct SphereV {
    float4 pos : SV_Position;
    nointerpolation float4 sph   : SPH;     // view-space centre, radius
    nointerpolation uint2  idv   : IDV;
    nointerpolation uint   state : STATE;
};

[shader("vertex")]
SphereV sphere_vs(uint vid : SV_VertexID, uint inst : SV_InstanceID) {
    Args     a = MD_ARGS;
    Batch    b = a.batches[a.slot];
    ViewData v = *a.view;
    SphereV  o;
    o.pos   = CULLED;
    o.sph   = float4(0.0);
    o.idv   = uint2(0u);
    o.state = 0u;

    uint e      = vid / 6u;
    uint corner = vid % 6u;
    if (e >= b.count) return o;
    uint pi = point_of(b, e, 0u, 1u);
    uint ai = attr_of(b, pi);
    if (!attr_visible(b, ai)) return o;
    uint st = attr_state(b, ai);
    if (a.mode == 1u && st == 0u) return o;

    Mat    T  = load_transform(b.transforms, inst);
    float  s  = length(T.c0.xyz);
    float  rad = sphere_radius(b, pi) * s;
    if (!(rad > 0.0)) return o;
    float3 cw = mulp(T, load_pos(b, pi));
    float3 c  = mulp(M(v, view), cw);

    // A quad on the plane tangent to the sphere at its point nearest the eye,
    // just large enough to hold the sphere's silhouette as seen through it.
    // The sphere lies wholly behind that plane, so every depth the fragment
    // stage writes is farther than the quad's: SV_DepthLessEqual holds.
    float3 dir;
    float  qr;
    if (v.ortho != 0u) {
        dir = float3(0.0, 0.0, -1.0);
        qr  = rad;
    } else {
        float d = length(c);
        if (d <= rad * 1.0001) return o;            // eye inside: not drawn
        dir = c / d;
        qr  = rad * (d - rad) / sqrt(d * d - rad * rad);
    }
    float3 np = c - dir * rad;
    float3 up = abs(dir.y) < 0.99 ? float3(0.0, 1.0, 0.0) : float3(1.0, 0.0, 0.0);
    float3 u  = normalize(cross(up, dir));
    float3 w  = cross(dir, u);
    float2 k;
    switch (corner) {
    case 0u: k = float2(-1.0, -1.0); break;
    case 1u: k = float2( 1.0, -1.0); break;
    case 2u: k = float2( 1.0,  1.0); break;
    case 3u: k = float2(-1.0, -1.0); break;
    case 4u: k = float2( 1.0,  1.0); break;
    default: k = float2(-1.0,  1.0); break;
    }
    float3 p = np + (u * k.x + w * k.y) * qr;

    o.pos   = mul4(M(v, proj), float4(p, 1.0));
    o.sph   = float4(c, rad);
    o.idv   = encode_id(a.slot, 0u, uint64_t(inst) * uint64_t(b.count) + uint64_t(e));
    o.state = st;
    return o;
}

struct VisOut {
    uint2 id    : SV_Target0;
    float depth : SV_DepthLessEqual;
};

struct MaskOut {
    float4 m    : SV_Target0;
    float depth : SV_DepthLessEqual;
};

[shader("fragment")]
VisOut sphere_fs(SphereV i) {
    ViewData v = *MD_ARGS.view;
    Ray   r = pixel_ray(v, i.pos.xy);
    float t = sphere_hit(r, i.sph.xyz, i.sph.w, false);
    if (t < 0.0) discard;
    VisOut o;
    o.id    = i.idv;
    o.depth = clip_depth(v, r.o + r.d * t);
    return o;
}

[shader("fragment")]
MaskOut sphere_mask_fs(SphereV i) {
    ViewData v = *MD_ARGS.view;
    Ray   r = pixel_ray(v, i.pos.xy);
    float t = sphere_hit(r, i.sph.xyz, i.sph.w, false);
    if (t < 0.0) discard;
    MaskOut o;
    o.m     = float4(0.0, (i.state & 1u) != 0u ? 1.0 : 0.0, 0.0, (i.state & 2u) != 0u ? 1.0 : 0.0);
    o.depth = clip_depth(v, r.o + r.d * t);
    return o;
}

struct CapsuleV {
    float4 pos : SV_Position;
    nointerpolation float4 pa    : PA;      // view space, w: radius
    nointerpolation float4 pb    : PB;
    nointerpolation uint2  idv   : IDV;
    nointerpolation uint   state : STATE;
};

// The box around a capsule, 12 triangles counter-clockwise seen from outside.
// Corner c is (c & 1, c >> 1 & 1, c >> 2 & 1) in the capsule's frame (u, v, axis).
static const uint BOX[36] = {
    0, 4, 6,  0, 6, 2,      // -u
    1, 3, 7,  1, 7, 5,      // +u
    0, 1, 5,  0, 5, 4,      // -v
    2, 6, 7,  2, 7, 3,      // +v
    0, 2, 3,  0, 3, 1,      // -axis
    4, 5, 7,  4, 7, 6,      // +axis
};

[shader("vertex")]
CapsuleV capsule_vs(uint vid : SV_VertexID, uint inst : SV_InstanceID) {
    Args     a = MD_ARGS;
    Batch    b = a.batches[a.slot];
    ViewData v = *a.view;
    CapsuleV o;
    o.pos   = CULLED;
    o.pa    = float4(0.0);
    o.pb    = float4(0.0);
    o.idv   = uint2(0u);
    o.state = 0u;

    uint e  = vid / 36u;
    uint vk = vid % 36u;
    if (e >= b.count) return o;
    uint p0 = point_of(b, e, 0u, 2u);
    uint p1 = point_of(b, e, 1u, 2u);
    uint a0 = attr_of(b, p0);
    uint a1 = attr_of(b, p1);
    if (!attr_visible(b, a0) || !attr_visible(b, a1)) return o;
    uint st = attr_state(b, a0) & attr_state(b, a1);
    if (a.mode == 1u && st == 0u) return o;

    Mat    T   = load_transform(b.transforms, inst);
    float  s   = length(T.c0.xyz);
    float  rad = capsule_radius(b, e) * s;
    if (!(rad > 0.0)) return o;
    float3 pa = mulp(M(v, view), mulp(T, load_pos(b, p0)));
    float3 pb = mulp(M(v, view), mulp(T, load_pos(b, p1)));

    float3 ax  = pb - pa;
    float  len = length(ax);
    float3 z   = len > 1e-20 ? ax / len : float3(0.0, 0.0, 1.0);
    float3 hlp = abs(z.x) < 0.9 ? float3(1.0, 0.0, 0.0) : float3(0.0, 1.0, 0.0);
    float3 x   = normalize(cross(z, hlp));
    float3 y   = cross(z, x);               // (x, y, z) right-handed

    // A mirroring transform would turn the box inside out: keep it outward.
    float det = dot(cross(T.c0.xyz, T.c1.xyz), T.c2.xyz);
    uint  tri = vk / 3u;
    uint  kk  = vk % 3u;
    if (det < 0.0 && kk != 0u) kk = 3u - kk;
    uint  c   = BOX[tri * 3u + kk];

    float3 p = pa
             + x * (((c & 1u) != 0u) ? rad : -rad)
             + y * (((c & 2u) != 0u) ? rad : -rad)
             + z * (((c & 4u) != 0u) ? len + rad : -rad);

    o.pos   = mul4(M(v, proj), float4(p, 1.0));
    o.pa    = float4(pa, rad);
    o.pb    = float4(pb, 0.0);
    o.idv   = encode_id(a.slot, 0u, uint64_t(inst) * uint64_t(b.count) + uint64_t(e));
    o.state = st;
    return o;
}

[shader("fragment")]
VisOut capsule_fs(CapsuleV i) {
    ViewData v = *MD_ARGS.view;
    Ray   r = pixel_ray(v, i.pos.xy);
    float t = capsule_hit(r, i.pa.xyz, i.pb.xyz, i.pa.w);
    if (t < 0.0) discard;
    VisOut o;
    o.id    = i.idv;
    o.depth = clip_depth(v, r.o + r.d * t);
    return o;
}

[shader("fragment")]
MaskOut capsule_mask_fs(CapsuleV i) {
    ViewData v = *MD_ARGS.view;
    Ray   r = pixel_ray(v, i.pos.xy);
    float t = capsule_hit(r, i.pa.xyz, i.pb.xyz, i.pa.w);
    if (t < 0.0) discard;
    MaskOut o;
    o.m     = float4(0.0, (i.state & 1u) != 0u ? 1.0 : 0.0, 0.0, (i.state & 2u) != 0u ? 1.0 : 0.0);
    o.depth = clip_depth(v, r.o + r.d * t);
    return o;
}

// ---------------------------------------------------------------------------
// Resolve
// ---------------------------------------------------------------------------

struct Hit {
    bool   valid;
    uint   batch;
    uint   instance;
    uint   element;
    float3 p_view;
    float3 n_view;
    float3 p_world;
    float3 p_world_prev;        // where that surface point was at the previous raster
    float4 color;               // sRGB-encoded, before tint
    uint   attr;                // of the nearest point
    uint   id;
};

float4 unpack_rgba8(uint c) {
    return float4(float(c & 0xFFu), float((c >> 8) & 0xFFu), float((c >> 16) & 0xFFu), float(c >> 24)) / 255.0;
}

float4 attr_color(Batch b, uint ai) {
    uint c;
    if (b.color_mode == COLOR_UNIFORM) {
        c = b.uniform_color;
    } else if (b.color_mode == COLOR_PALETTE) {
        uint idx = 0u;
        if ((b.flags & F_PALETTE_IDX) != 0u) {
            uint w = b.palette_index[ai >> 1];
            idx = ((ai & 1u) != 0u) ? (w >> 16) : (w & 0xFFFFu);
        }
        if (idx >= b.palette_count) idx = 0u;
        c = b.palette_count > 0u ? b.palette[idx] : 0xFFFFFFFFu;
    } else {
        c = ((b.flags & F_COLOR) != 0u) ? b.color[ai] : 0xFFFFFFFFu;
    }
    return unpack_rgba8(c);
}

float3 srgb_to_linear(float3 c) {
    return select(c <= 0.04045, c / 12.92, pow((c + 0.055) / 1.055, 2.4));
}

float3 linear_to_srgb(float3 c) {
    return select(c <= 0.0031308, c * 12.92, 1.055 * pow(c, 1.0 / 2.4) - 0.055);
}

float2 oct_encode(float3 n) {
    n /= (abs(n.x) + abs(n.y) + abs(n.z));
    float2 p = n.xy;
    if (n.z < 0.0) {
        float2 s = select(p >= 0.0, float2(1.0), float2(-1.0));
        p = (1.0 - abs(p.yx)) * s;
    }
    return p;
}

Hit eval_hit(ViewData v, Batch* batches, uint2 idv, float2 pix) {
    Hit h;
    h.valid = false;
    uint slot1 = idv.y >> 20;
    if (slot1 == 0u) return h;
    h.batch = slot1 - 1u;
    Batch b = batches[h.batch];
    uint64_t linear = (uint64_t(idv.y & 0x3FFFFu) << 32) | uint64_t(idv.x);
    uint cnt = max(b.count, 1u);
    h.instance = uint(linear / uint64_t(cnt));
    h.element  = uint(linear - uint64_t(h.instance) * uint64_t(cnt));
    h.id       = ((b.flags & F_ID) != 0u) ? b.id[h.element] : h.element;
    h.attr     = 0xFFFFFFFFu;
    if (b.kind == KIND_EXTERNAL) return h;

    Mat   T  = load_transform(b.transforms, h.instance);
    Mat   Tp = load_transform(b.transforms_prev, h.instance);
    float s  = length(T.c0.xyz);
    Ray   r  = pixel_ray(v, pix);

    if (b.kind == KIND_SPHERES) {
        uint   pi  = point_of(b, h.element, 0u, 1u);
        float3 lc  = load_pos(b, pi);
        float3 cw  = mulp(T, lc);
        float3 cwp = mulp(Tp, lc - load_vel(b, pi) * v.dt);
        float3 c   = mulp(M(v, view), cw);
        float  rad = sphere_radius(b, pi) * s;
        float  t   = sphere_hit(r, c, rad, true);
        h.p_view   = r.o + r.d * t;
        h.n_view   = normalize(h.p_view - c);
        h.p_world  = mulp(M(v, inv_view), h.p_view);
        h.p_world_prev = h.p_world - (cw - cwp);
        h.attr     = attr_of(b, pi);
        h.color    = attr_color(b, h.attr);
        h.valid    = true;
    } else if (b.kind == KIND_CAPSULES) {
        uint   p0  = point_of(b, h.element, 0u, 2u);
        uint   p1  = point_of(b, h.element, 1u, 2u);
        float3 l0  = load_pos(b, p0);
        float3 l1  = load_pos(b, p1);
        float3 w0  = mulp(T, l0);
        float3 w1  = mulp(T, l1);
        float3 w0p = mulp(Tp, l0 - load_vel(b, p0) * v.dt);
        float3 w1p = mulp(Tp, l1 - load_vel(b, p1) * v.dt);
        float3 pa  = mulp(M(v, view), w0);
        float3 pb  = mulp(M(v, view), w1);
        float  rad = capsule_radius(b, h.element) * s;
        float  t   = capsule_hit(r, pa, pb, rad);
        if (t < 0.0) {
            // A rounding miss at the silhouette: the nearest point of the ray
            // to the axis stands in.
            float3 ba = pb - pa;
            float3 oa = r.o - pa;
            float  bd = dot(ba, r.d);
            float  den = max(dot(ba, ba) - bd * bd, 1e-20);
            t = max((dot(oa, ba) * bd - dot(oa, r.d) * dot(ba, ba)) / den, 0.0);
        }
        h.p_view   = r.o + r.d * t;
        float k    = capsule_param(h.p_view, pa, pb);
        h.n_view   = normalize(h.p_view - lerp(pa, pb, k));
        h.p_world  = mulp(M(v, inv_view), h.p_view);
        h.p_world_prev = h.p_world - lerp(w0 - w0p, w1 - w1p, k);
        uint   a0  = attr_of(b, p0);
        uint   a1  = attr_of(b, p1);
        h.attr     = k < 0.5 ? a0 : a1;
        float  sh  = saturate(b.sharpness);
        float  kk  = saturate((k - 0.5) / max(1.0 - sh, 1e-3) + 0.5);
        h.color    = lerp(attr_color(b, a0), attr_color(b, a1), kk);
        h.valid    = true;
    }
    return h;
}

float2 to_pixels(ViewData v, float4 clip) {
    float2 ndc = clip.xy / clip.w;
    return float2((ndc.x * 0.5 + 0.5) * v.viewport.x, (0.5 - ndc.y * 0.5) * v.viewport.y);
}

[shader("compute")]
[numthreads(8, 8, 1)]
void resolve(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    if (tid.x >= a.width || tid.y >= a.height) return;
    int2  px  = int2(tid.xy);
    uint2 idv = a.id.Load(int3(px, 0));
    if ((idv.y >> 20) == 0u) {
        a.color[px]  = float4(0.0);
        a.normal[px] = float2(0.0);
        a.motion[px] = float2(0.0);
        return;
    }
    ViewData v = *a.view;
    Hit h = eval_hit(v, a.batches, idv, float2(tid.xy) + 0.5);
    if (!h.valid) return;       // EXTERNAL: the caller's

    Batch  b   = a.batches[h.batch];
    float3 lin = srgb_to_linear(h.color.rgb);
    float  g   = dot(lin, float3(0.2126, 0.7152, 0.0722));
    lin = lerp(lin, float3(g), saturate(b.desaturate));
    lin = lerp(lin, lin * b.tint.rgb, saturate(b.tint.a));
    a.color[px]  = float4(linear_to_srgb(saturate(lin)), 1.0);
    a.normal[px] = oct_encode(h.n_view);

    float4 cur  = mul4(M(v, proj_cur),  mul4(M(v, view),      float4(h.p_world, 1.0)));
    float4 prev = mul4(M(v, proj_prev), mul4(M(v, view_prev), float4(h.p_world_prev, 1.0)));
    a.motion[px] = to_pixels(v, cur) - to_pixels(v, prev);
}

// The visible part of the masks, from the ids: r state 0, b state 1. The
// raster pass that follows adds g and a (anywhere, hidden too) by max.
[shader("compute")]
[numthreads(8, 8, 1)]
void mask_visible(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    if (tid.x >= a.width || tid.y >= a.height) return;
    int2  px  = int2(tid.xy);
    uint2 idv = a.id.Load(int3(px, 0));
    float4 m  = float4(0.0);
    uint slot1 = idv.y >> 20;
    if (slot1 != 0u) {
        Batch b = a.batches[slot1 - 1u];
        if (b.kind == KIND_SPHERES || b.kind == KIND_CAPSULES) {
            uint64_t linear = (uint64_t(idv.y & 0x3FFFFu) << 32) | uint64_t(idv.x);
            uint cnt = max(b.count, 1u);
            uint e   = uint(linear % uint64_t(cnt));
            uint st;
            if (b.kind == KIND_SPHERES) {
                st = attr_state(b, attr_of(b, point_of(b, e, 0u, 1u)));
            } else {
                st = attr_state(b, attr_of(b, point_of(b, e, 0u, 2u))) & attr_state(b, attr_of(b, point_of(b, e, 1u, 2u)));
            }
            m = float4((st & 1u) != 0u ? 1.0 : 0.0, 0.0, (st & 2u) != 0u ? 1.0 : 0.0, 0.0);
        }
    }
    a.mask[px] = m;
}

// One pixel, decoded: id.x, id.y, element id, attr, position xyz, depth.
[shader("compute")]
[numthreads(1, 1, 1)]
void pick(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    uint* o = a.pick_out;
    for (uint i = 0u; i < 12u; ++i) o[i] = 0u;
    if (a.px >= a.width || a.py >= a.height) return;
    uint2 idv = a.id.Load(int3(int(a.px), int(a.py), 0));
    o[0] = idv.x;
    o[1] = idv.y;
    if ((idv.y >> 20) == 0u) return;
    ViewData v = *a.view;
    Hit h = eval_hit(v, a.batches, idv, float2(float(a.px), float(a.py)) + 0.5);
    o[2] = h.id;
    o[3] = h.attr;
    if (h.valid) {
        o[4] = asuint(h.p_world.x);
        o[5] = asuint(h.p_world.y);
        o[6] = asuint(h.p_world.z);
        o[7] = asuint(clip_depth(v, h.p_view));
    }
}
