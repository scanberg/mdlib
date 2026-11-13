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
    uint   total_work;          // elements times instances over all batches
    float  near_z;              // view-space z of the near plane (negative)
    uint   hiz_valid;           // the pyramid holds the previous raster at this size
    uint   hiz_w;               // its level 0 (half the target, rounded up)
    uint   hiz_h;
    uint   hiz_levels;
    uint   occlusion;           // cull against the pyramid
    uint   compute_raster;      // rasterise small primitives in compute
    float  small_px;            // up to this many pixels across counts as small
    uint   total_clusters;      // clusters over all batches (refit)
    uint   total_cluster_work;  // clusters times instances (cluster cull)
    uint   refit_mode;          // 0 bounds and state, 1 state only
    float  jitter_px;           // how far this raster's jitter moves the image (pixels)
    float  occl_margin;         // phase 1: occluded only when nearer by this many radii too
    float  occl_min_px;         // elements no wider are not tested for occlusion
    uint   _vpad1;
    uint   _vpad2;
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
static const uint F_CLUSTERS     = 1u << 12;  // explicit cluster offsets
static const uint F_SHARED       = 1u << 13;  // bounds shared with an earlier batch: no refit
static const uint F_RADIUS_TABLE = 1u << 14;  // `radius` holds a type per point (low 16 bits; high: flags)
static const uint F_FADE_OUT     = 1u << 15;  // fades out across a far band: fade_out_a, fade_out_b
static const uint F_FADE_IN      = 1u << 16;  // fades in across a near band: fade_in_a, fade_in_b
static const uint F_FADE_SHRINK  = 1u << 17;  // ... its sizes scaled
static const uint F_FADE_THIN    = 1u << 18;  // ... its elements dropped one by one
static const uint F_FADE         = F_FADE_OUT | F_FADE_IN;
static const uint F_SPL_AXIS     = 1u << 19;  // SPLINES: spl_axis
static const uint F_SPL_PROFILE  = 1u << 20;  // SPLINES: spl_profile
static const uint F_SPL_THROUGH  = 1u << 21;  // SPLINES: through their control points (spline_controls)

// md_render_kind_t
static const uint KIND_SPHERES  = 1u;
static const uint KIND_CAPSULES = 2u;
static const uint KIND_POLYGONS = 3u;
static const uint KIND_SPLINES  = 4u;
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
    uint    work_offset;        // first work item (element x instance) of the state cull
    uint    rad_stride;         // in bytes

    // Clusters: the unit of the first level of the cull. Cluster c holds the
    // elements cluster_offset[c] .. cluster_offset[c + 1], or 64c .. 64c + 64
    // without explicit offsets. Bounds are in the batch's own space, of the
    // element centres (spheres) or ends (capsules); each batch pads them by
    // its own style at cull time, so batches over the same points share them.
    uint*   cluster_offset;
    float4* bounds;             // [2 per cluster]: min.xyz, max radius factor; max.xyz, unused
    uint*   cstate;             // [2 per cluster]: OR of the state bits of its elements, visible elements
    float*  radius_table;       // F_RADIUS_TABLE: float per type
    uint    cluster_count;
    uint    cluster_work_offset;// first cluster work item (cluster x instance)
    uint    cluster_base;       // first cluster over all batches (refit)
    uint    radius_table_count;

    // The fade bands for this camera (md_render.c, fade_bands): how far
    // through a band an element at view-space depth d is, saturate(a d + b)
    float   fade_out_a;
    float   fade_out_b;
    float   fade_in_a;
    float   fade_in_b;

    float4* spl_axis;           // SPLINES: per control point, the profile's long axis (xyz)
    float4* spl_profile;        // SPLINES: per control point, its part's half widths and radii at its start and end
    uint    poly_n;             // POLYGONS: corners of each
    uint    _pad_b0;            // to the float4 alignment
    uint    _pad_b1;
    uint    _pad_b2;
};

// The per-view lists the cull writes and everything after it reads. Cluster
// entries: the batch's local cluster work index and (slot + 1) << 20. Element
// entries: the element's global work index (its batch's work_offset plus
// instance * count + element; decode_work), half the size of an id.
struct Lists {
    uint2*    cvis;             // clusters to expand: phase 1, then phase 2 after them
    uint2*    cocc;             // clusters occluded in phase 1
    uint2*    cst;              // clusters with a state bit, in the frustum
    uint*     hw0;              // spheres for the hardware raster
    uint*     hw1;              // capsules
    uint*     st0;              // spheres with a state bit, for the masks
    uint*     st1;              // capsules with a state bit
    uint*     small;            // primitives for the compute raster
    uint*     occ;              // occluded in phase 1, tested again in phase 2
    uint*     counters;         // see C_* below
    uint64_t* vis64;            // compute raster: depth << 32 | index into small
    uint*     hw2;              // spline parts for the hardware raster (they have no compute raster)
    uint*     st2;              // spline parts with a state bit
    uint*     hw3;              // polygons for the hardware raster (nor have they)
    uint*     st3;              // polygons with a state bit
};

// counters[]: counts, then the indirect commands made of them
static const uint C_HW1        = 0;    // + kind index (0 spheres, 1 capsules)
static const uint C_HW2        = 2;
static const uint C_ST         = 4;
static const uint C_SMALL1     = 6;
static const uint C_SMALL2     = 7;
static const uint C_OCC        = 8;
static const uint C_HIDDEN     = 9;
static const uint C_OUTSIDE    = 10;
static const uint C_OCCLUDED   = 11;
static const uint C_CVIS1      = 12;
static const uint C_CVIS2      = 13;
static const uint C_COCC       = 14;
static const uint C_CST        = 15;
static const uint C_GRID_CVIS  = 52;   // [phase - 1], 4 words apart: 52, 56
static const uint C_GRID_COCC  = 60;
static const uint C_GRID_CST   = 64;
static const uint C_GRID_STATE = 68;   // unused
static const uint C_CL_VIS     = 72;   // clusters drawn from, phase 1 and 2: 72, 73
static const uint C_CL_OCC     = 74;   // clusters occluded after both phases
static const uint C_CMD        = 16;   // md_gpu_draw_cmd_t [phase - 1][kind]: 16..31
static const uint C_STCMD      = 32;   // [kind]: 32..39
static const uint C_GRID_SMALL = 40;   // md_gpu_grid_t [phase - 1], 4 words apart: 40, 44
static const uint C_GRID_OCC   = 48;
static const uint C_ICMD       = 80;   // md_gpu_draw_indexed_cmd_t, capsules [phase - 1]: 80..89
static const uint C_STICMD     = 90;   // the state list's capsules: 90..94
static const uint C_LOD        = 95;   // elements dropped by their batch's fade
static const uint C_HWS1       = 96;   // spline parts for the hardware raster, phase 1
static const uint C_HWS2       = 97;   // and phase 2
static const uint C_STS        = 98;   // spline parts with a state bit
static const uint C_ICMDS      = 100;  // md_gpu_draw_indexed_cmd_t, spline parts [phase - 1]: 100..109
static const uint C_STICMDS    = 110;  // the state list's spline parts: 110..114
static const uint C_HWP1       = 116;  // polygons for the hardware raster, phase 1
static const uint C_HWP2       = 117;  // and phase 2
static const uint C_STP        = 118;  // polygons with a state bit
static const uint C_ICMDP      = 120;  // md_gpu_draw_indexed_cmd_t, polygons [phase - 1]: 120..129
static const uint C_STICMDP    = 130;  // the state list's polygons: 130..134
static const uint C_COUNT      = 136;

// One argument struct for every entry point (md_gpu shares one root between a
// draw's two stages, and keeping kernels on it too means one C mirror).
//   raster:      view, batches, lists, mode (0 visibility, 1 mask), list (0..7:
//                hw0, hw1, st0, st1, hw2, st2, hw3, st3)
//   cull, finalize, raster_small: view, batches, lists, mode = phase (1, 2)
//   resolve, mask_visible, pick: view, batches, the targets, size, px, py, pick_out
//   hiz_reduce:  src (mip src_mip, src_w x src_h), dst (dst_w x dst_h)
struct Args {
    ViewData* view;
    Batch*    batches;
    uint      mode;
    uint      list;
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
    Lists*    lists;
    DescriptorHandle<Texture2D<float>>     src;
    DescriptorHandle<Texture2D<float>>     hiz;
    DescriptorHandle<RWTexture2D<float>>   dst;
    uint      src_mip;
    uint      src_w;
    uint      src_h;
    uint      dst_w;
    uint      dst_h;
    uint      _pad0;
};
MD_SHADER_ARGS(Args);

// ---------------------------------------------------------------------------
// Fetch
// ---------------------------------------------------------------------------

// No pointer of a batch is null (md_render.c, no_null): an absent array
// points at a block of zeros, and its index is gated by its flag to 0. A load
// executed in a branch not taken reads that block. (lavapipe executes such
// loads at its first lane's address whatever the mask; a driver that
// speculates loads out of branches would too.)
uint gated(Batch b, uint flag, uint i) { return (b.flags & flag) != 0u ? i : 0u; }

float3 load_pos(Batch b, uint i) {
    uint* p = b.position + i * (b.stride >> 2);
    return float3(asfloat(p[0]), asfloat(p[1]), asfloat(p[2]));
}

float3 load_vel(Batch b, uint i) {
    if ((b.flags & F_VELOCITY) == 0u) return float3(0.0);
    uint* p = b.velocity + gated(b, F_VELOCITY, i) * (b.vel_stride >> 2);
    uint w0 = p[0];
    uint w1 = p[1];
    return float3(f16tof32(w0 & 0xFFFFu), f16tof32(w0 >> 16), f16tof32(w1 & 0xFFFFu));
}

uint attr_of(Batch b, uint pi) {
    return ((b.flags & F_ATTR_INDEX) != 0u) ? b.attr_index[gated(b, F_ATTR_INDEX, pi)] : pi;
}

bool test_bit(uint* set, uint i) {
    return ((set[i >> 5] >> (i & 31u)) & 1u) != 0u;
}

bool attr_visible(Batch b, uint a) {
    if ((b.flags & F_VISIBLE0) != 0u && !test_bit(b.visible0, gated(b, F_VISIBLE0, a))) return false;
    if ((b.flags & F_VISIBLE1) != 0u && !test_bit(b.visible1, gated(b, F_VISIBLE1, a))) return false;
    return true;
}

uint attr_state(Batch b, uint a) {
    uint s = 0u;
    if ((b.flags & F_STATE0) != 0u && test_bit(b.state0, gated(b, F_STATE0, a))) s |= 1u;
    if ((b.flags & F_STATE1) != 0u && test_bit(b.state1, gated(b, F_STATE1, a))) s |= 2u;
    return s;
}

// Point k (0 or 1) of element e with `per` points per element
uint point_of(Batch b, uint e, uint k, uint per) {
    uint i = e * per + k;
    return ((b.flags & F_INDEX) != 0u) ? b.index[gated(b, F_INDEX, i)] : i;
}

Mat load_transform(float4* t, uint inst) {
    Mat m;
    m.c0 = t[inst * 4u + 0u];
    m.c1 = t[inst * 4u + 1u];
    m.c2 = t[inst * 4u + 2u];
    m.c3 = t[inst * 4u + 3u];
    return m;
}

// A point's radius: its own, or its type's (F_RADIUS_TABLE; the type comes in
// the same load as the position when the points are float4). The lookup read
// from memory costs the GTX 1060's refit at 10M atoms 0.04 ms over a float
// radius (md_render_bench, "-- radii"). Copied into group memory first it
// cost 0.1 ms more (the copy and a barrier per group), and the kernels that
// did it took more registers (the refit 56 -> 64, the cull and the compute
// raster 80 -> 96).
float point_radius(Batch b, uint pi) {
    if ((b.flags & F_RADIUS) == 0u) return 0.0;
    float v = b.radius[gated(b, F_RADIUS, pi) * (b.rad_stride >> 2)];
    if ((b.flags & F_RADIUS_TABLE) == 0u) return v;
    uint t = asuint(v) & 0xFFFFu;
    return t < b.radius_table_count ? b.radius_table[t < b.radius_table_count ? t : 0u] : 0.0;
}

float sphere_radius(Batch b, uint pi) {
    float r = point_radius(b, pi);
    return r * b.s0 + b.s1;
}

float capsule_radius(Batch b, uint e) {
    float r = b.s0;
    if ((b.flags & F_RADIUS_SCALE) != 0u) r *= b.radius_scale[gated(b, F_RADIUS_SCALE, e)];
    return r;
}

// ---------------------------------------------------------------------------
// Rays and intersections (view space)
// ---------------------------------------------------------------------------

struct Ray { float3 o; float3 d; };

// Every multiply-add explicit: the compute raster and the fragment shaders
// must make the same ray to the bit, and under fast math a compiler may fuse
// a * b + c in one stage and not in another (a jittered projection's inverse
// has terms that make the two differ)
float3 unproject(Mat inv_proj, float2 ndc, float z) {
    float4 p = fma(inv_proj.c0, float4(ndc.x), fma(inv_proj.c1, float4(ndc.y), fma(inv_proj.c2, float4(z), inv_proj.c3)));
    return p.xyz / p.w;
}

// The ray of pixel position `pix` (pixel units, top-left origin; a pixel's
// centre is at +0.5), through the projection as it was rasterised.
Ray pixel_ray(ViewData v, float2 pix) {
    float2 ndc = float2(fma(pix.x, v.viewport.z * 2.0, -1.0), fma(-pix.y, v.viewport.w * 2.0, 1.0));
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

// Where a sphere and a capsule are the same surface (the joints of licorice,
// whose spheres have the sticks' radius), the sphere must win whatever path
// rasterised each and in whatever order. Capsule depths are one ulp farther
// (reverse Z: smaller), which breaks exactly those ties and nothing visible.
float capsule_depth(float d) {
    return d > 0.0 ? asfloat(asuint(d) - 1u) : d;
}

// ---------------------------------------------------------------------------
// Ids
// ---------------------------------------------------------------------------

uint2 encode_id(uint slot, uint part, uint64_t linear) {
    uint hi = uint(linear >> 32) & 0x3FFFFu;
    return uint2(uint(linear & 0xFFFFFFFFull), ((slot + 1u) << 20) | ((part & 3u) << 18) | hi);
}

// ---------------------------------------------------------------------------
// Elements
// ---------------------------------------------------------------------------

struct Elem {
    uint slot;
    uint inst;
    uint e;
};

Elem decode(Batch* batches, uint2 ent) {
    Elem el;
    el.slot = (ent.y >> 20) - 1u;
    uint64_t linear = (uint64_t(ent.y & 0x3FFFFu) << 32) | uint64_t(ent.x);
    uint cnt = max(batches[el.slot].count, 1u);
    el.inst = uint(linear / uint64_t(cnt));
    el.e    = uint(linear - uint64_t(el.inst) * uint64_t(cnt));
    return el;
}

// An element in view space: spheres a = b = centre; capsules a, b the ends.
// Splines and polygons are heavy code: in a kernel that also takes spheres
// and capsules they would set its registers, and so its occupancy, for all
// of them (splines took resolve from 168 to 242 on Pascal, the compute
// raster from 80 to 168). So the kernels that run per element or per pixel
// come in twins: the main one takes spheres and capsules, its twin (SPL 1)
// splines and polygons, and the twin is dispatched only for frames that
// have those. Both walk the same work; each passes over what is the other's.
bool heavy_kind(uint kind) { return kind == KIND_SPLINES || kind == KIND_POLYGONS; }
bool kind_taken<let SPL : int>(uint kind) {
    return SPL != 0 ? heavy_kind(kind) : !heavy_kind(kind);
}

struct Geo {
    bool   ok;                  // visible, of positive radius, not faded
    bool   faded;               // would be drawn but for its batch's fade
    uint   kind;                // the hardware raster's: 0 spheres, 1 capsules, 2 spline parts, 3 polygons
    uint   st;                  // state bits
    float3 a;
    float3 b;
    float  r;
};

// Fading with distance (md_render_batch_t.fade). A band is a ramp in
// view-space depth: t = saturate(a d + b), 0 where it begins (nearer), 1
// where it ends. Across the far band a batch fades out, across the near one
// in. The host gives the two sides of a handover the same a and b, and each
// side computes t by the same expression: they see the same t to the bit.
float fade_t(float a, float b, float d) { return saturate(fma(a, d, b)); }

// How far through its far band (0 without one) and its near band (1 without
// one) an element at depth d is
float2 fade_ts(Batch b, float d) {
    float to = (b.flags & F_FADE_OUT) != 0u ? fade_t(b.fade_out_a, b.fade_out_b, d) : 0.0;
    float ti = (b.flags & F_FADE_IN)  != 0u ? fade_t(b.fade_in_a,  b.fade_in_b,  d) : 1.0;
    return float2(to, ti);
}

// SHRINK: the size factor
float fade_size(float2 t) { return (1.0 - t.x) * t.y; }

// THIN: an attribute is kept while its hash h (in [0, 1), 24 bits) is at
// least how far through the far band and below how far through the near one.
// Across a shared band the side fading out keeps h >= t and the side fading
// in h < t: each attribute on exactly one side.
uint fade_hash(uint x) {
    x ^= x >> 16; x *= 0x7FEB352Du;
    x ^= x >> 15; x *= 0x846CA68Bu;
    x ^= x >> 16;
    return x;
}
bool fade_keeps(float2 t, uint attr) {
    float h = float(fade_hash(attr) >> 8) * (1.0 / 16777216.0);
    return h >= t.x && h < t.y;
}

// ---------------------------------------------------------------------------
// Splines (md_render_spline_style_t): a stadium swept along a uniform cubic
// B-spline. Element e is control point e and draws its part: the curve
// parameter u in [-0.5, 0.5] around it (u = 0 at e), [0, 0.5] for the first
// of a spline and [-0.5, 0] for the last. The part is cut into K pieces
// (more the more it bends on screen); a piece is the ball of the radius
// around the segments between its two cross sections: the rims, two
// capsules along the ends of the cross sections, and the faces, two
// bilinear patches at the radius on either side. Neighbouring parts compute their shared
// cross section from the same control points with the same expressions:
// they meet exactly. Every stage (cull, both rasters, resolve, pick, mask)
// builds the part with spline_part and cuts it with spline_pieces, so they
// agree to the bit.
// ---------------------------------------------------------------------------

static const uint SPLINE_FIRST = 0x80000000u;
static const uint SPLINE_LAST  = 0x40000000u;
static const uint SPLINE_POINT = 0x3FFFFFFFu;

struct SplinePart {
    bool   ok;
    float3 p0, p1, p2, p3, p4;  // control points e-2 .. e+2 in view space, mirrored past the ends
    float3 w0, w1, w2, w3, w4;  // their axes, unit, the sign the caller gave
    float  u0, u1;              // the part's range of u
    float  a0, a1;              // half widths at u0 and u1 (scaled, shrunk)
    float  r0, r1;              // radii at u0 and u1 (scaled, shrunk)
    float  ru0, ru1;            // the radius's rate along u at u0 and u1: at
                                // a joint without a cap the mean of the two
                                // parts' (spline_frame's tilt, continuous)
    uint   caps;                // 1: a cap at u0, 2: at u1
};

float3 spline_axis_of(Batch b, uint e) {
    return (b.flags & F_SPL_AXIS) != 0u ? b.spl_axis[gated(b, F_SPL_AXIS, e)].xyz : float3(0.0);
}

// e's axis in view space, of unit length (0 stays 0): the blend weighs every
// control point alike, whatever the length the caller gave
float3 spline_axis_view(Mat V, Mat T, Batch b, uint e) {
    float3 w  = mulv(V, mulv(T, spline_axis_of(b, e)));
    float  ww = dot(w, w);
    return ww > 1e-30 ? w * rsqrt(ww) : float3(0.0);
}

// e's part's profile: half width at its start and end, radius at its start
// and end; without a profile the batch's flat (s1) and radius (s0)
float4 spline_profile_of(Batch b, uint e) {
    return (b.flags & F_SPL_PROFILE) != 0u ? b.spl_profile[gated(b, F_SPL_PROFILE, e)] : float4(b.s1, b.s1, b.s0, b.s0);
}

bool spline_visible(Batch b, uint ie) {
    return attr_visible(b, attr_of(b, ie & SPLINE_POINT));
}

// The control points of e's part, batch space: e - 2 .. e + 2, mirrored past
// the spline's ends (where the B-spline then ends at its end point). Through
// its control points (F_SPL_THROUGH) each is moved by one step of the
// inverse of the B-spline's smoothing, (q- + 4 q + q+) / 6: c = (8 q - q- -
// q+) / 6, which leaves an offset of the fourth difference over 36 instead
// of the second over 6 (md_gl's tubes.through). An end
// stays where it is (its mirrored neighbour makes the step nothing). Every
// load past an end behind an if and its index gated (no_null in md_render.c).
// The caller has checked that e is not alone.
struct SplineCtl {
    float3 p0, p1, p2, p3, p4;
    bool   hp, hn, hpp, hnn;            // e - 1, e + 1, e - 2, e + 2 are of the spline
    uint   im, ip;                      // the index entries of e - 1, e + 1 (0 past an end)
};

SplineCtl spline_controls(Batch b, uint e, uint ie) {
    SplineCtl k;
    k.hp = (ie & SPLINE_FIRST) == 0u;
    k.hn = (ie & SPLINE_LAST) == 0u;
    k.im = 0u; k.ip = 0u;
    if (k.hp) k.im = b.index[k.hp ? e - 1u : e];
    if (k.hn) k.ip = b.index[k.hn ? e + 1u : e];
    k.hpp = k.hp && (k.im & SPLINE_FIRST) == 0u;
    k.hnn = k.hn && (k.ip & SPLINE_LAST) == 0u;
    uint imm = 0u, ipp = 0u;
    if (k.hpp) imm = b.index[k.hpp ? e - 2u : e];
    if (k.hnn) ipp = b.index[k.hnn ? e + 2u : e];
    k.p2 = load_pos(b, ie & SPLINE_POINT);
    k.p0 = k.p1 = k.p3 = k.p4 = float3(0.0);
    if (k.hp) k.p1 = load_pos(b, k.im & SPLINE_POINT);
    if (k.hn) k.p3 = load_pos(b, k.ip & SPLINE_POINT);
    if (!k.hp) k.p1 = 2.0 * k.p2 - k.p3;
    if (!k.hn) k.p3 = 2.0 * k.p2 - k.p1;
    if (k.hpp) k.p0 = load_pos(b, imm & SPLINE_POINT);
    else       k.p0 = 2.0 * k.p1 - k.p2;
    if (k.hnn) k.p4 = load_pos(b, ipp & SPLINE_POINT);
    else       k.p4 = 2.0 * k.p3 - k.p2;
    if ((b.flags & F_SPL_THROUGH) != 0u) {
        bool   hppp = k.hpp && (imm & SPLINE_FIRST) == 0u;
        bool   hnnn = k.hnn && (ipp & SPLINE_LAST) == 0u;
        float3 m3 = float3(0.0), n3 = float3(0.0);
        if (hppp) m3 = load_pos(b, b.index[hppp ? e - 3u : e] & SPLINE_POINT);
        if (hnnn) n3 = load_pos(b, b.index[hnnn ? e + 3u : e] & SPLINE_POINT);
        float3 c2 = (k.hp && k.hn)  ? (8.0 * k.p2 - k.p1 - k.p3) / 6.0 : k.p2;
        float3 c1 = (k.hp && k.hpp) ? (8.0 * k.p1 - k.p0 - k.p2) / 6.0 : k.p1;
        float3 c3 = (k.hn && k.hnn) ? (8.0 * k.p3 - k.p2 - k.p4) / 6.0 : k.p3;
        float3 c0 = (k.hpp && hppp) ? (8.0 * k.p0 - m3 - k.p1) / 6.0 : k.p0;
        float3 c4 = (k.hnn && hnnn) ? (8.0 * k.p4 - k.p3 - n3) / 6.0 : k.p4;
        if (!k.hp)  c1 = 2.0 * c2 - c3;     // mirrored again, of the moved ones
        if (!k.hn)  c3 = 2.0 * c2 - c1;
        if (!k.hpp) c0 = 2.0 * c1 - c2;
        if (!k.hnn) c4 = 2.0 * c3 - c2;
        k.p0 = c0; k.p1 = c1; k.p2 = c2; k.p3 = c3; k.p4 = c4;
    }
    return k;
}

// Control point e's part, in view space. fade: the batch's SHRINK applies
// (by e's depth, as element_geo's test).
SplinePart spline_part(ViewData v, Batch b, uint inst, uint e, bool fade) {
    SplinePart sp;
    sp.ok = false;
    sp.p0 = sp.p1 = sp.p2 = sp.p3 = sp.p4 = float3(0.0);
    sp.w0 = sp.w1 = sp.w2 = sp.w3 = sp.w4 = float3(0.0);
    sp.u0 = sp.u1 = 0.0; sp.a0 = sp.a1 = 0.0; sp.r0 = sp.r1 = 0.0; sp.ru0 = sp.ru1 = 0.0; sp.caps = 0u;
    uint ie = b.index[e];
    bool hp = (ie & SPLINE_FIRST) == 0u;
    bool hn = (ie & SPLINE_LAST) == 0u;
    if (!hp && !hn) return sp;                          // alone: no curve
    SplineCtl k = spline_controls(b, e, ie);
    uint im = k.im, ip = k.ip;
    bool hpp = k.hpp, hnn = k.hnn;
    Mat  T = load_transform(b.transforms, inst);
    Mat  V = M(v, view);
    float s = length(T.c0.xyz);
    sp.p0 = mulp(V, mulp(T, k.p0));
    sp.p1 = mulp(V, mulp(T, k.p1));
    sp.p2 = mulp(V, mulp(T, k.p2));
    sp.p3 = mulp(V, mulp(T, k.p3));
    sp.p4 = mulp(V, mulp(T, k.p4));

    sp.w2 = spline_axis_view(V, T, b, e);
    sp.w1 = sp.w2; sp.w3 = sp.w2;
    if (hp) sp.w1 = spline_axis_view(V, T, b, e - 1u);
    if (hn) sp.w3 = spline_axis_view(V, T, b, e + 1u);
    sp.w0 = sp.w1; sp.w4 = sp.w3;
    if (hpp) sp.w0 = spline_axis_view(V, T, b, e - 2u);
    if (hnn) sp.w4 = spline_axis_view(V, T, b, e + 2u);

    float4 pf = spline_profile_of(b, e);
    float  f  = 1.0;
    if (fade && (b.flags & F_FADE) != 0u && (b.flags & F_FADE_SHRINK) != 0u) f = fade_size(fade_ts(b, -sp.p2.z));
    sp.u0 = hp ? -0.5 : 0.0;
    sp.u1 = hn ?  0.5 : 0.0;
    sp.a0 = max(pf.x, 0.0) * s * f;
    sp.a1 = max(pf.y, 0.0) * s * f;
    sp.r0 = max(pf.z, 0.0) * s * f;
    sp.r1 = max(pf.w, 0.0) * s * f;
    // caps: at the ends of the spline, next to a hidden part, at a step of
    // the profile. Where the joint is smooth, the radius's rate there is the
    // mean of the two parts' (a part's u spans 1, or 0.5 at a spline's end):
    // each part tilts its faces by it (spline_frame), and with a rate of its
    // own each would end its faces at its own tangent point on the joint's
    // balls, a sliver apart
    float rr = (pf.w - pf.z) / (sp.u1 - sp.u0);
    sp.ru0 = rr;
    sp.ru1 = rr;
    if (!hp) sp.caps |= 1u;
    else {
        float4 pm = spline_profile_of(b, e - 1u);
        if (!spline_visible(b, im) || any(pm.yw != pf.xz)) sp.caps |= 1u;
        else sp.ru0 = 0.5 * (rr + (pm.w - pm.z) / ((im & SPLINE_FIRST) != 0u ? 0.5 : 1.0));
    }
    if (!hn) sp.caps |= 2u;
    else {
        float4 pn = spline_profile_of(b, e + 1u);
        if (!spline_visible(b, ip) || any(pn.xz != pf.yw)) sp.caps |= 2u;
        else sp.ru1 = 0.5 * (rr + (pn.w - pn.z) / ((ip & SPLINE_LAST) != 0u ? 0.5 : 1.0));
    }
    sp.ru0 *= s * f;
    sp.ru1 *= s * f;
    sp.ok = max(sp.r0, sp.r1) > 0.0;
    return sp;
}

// The cross section at u: centre, tangent (unnormalised) and unit axis
// orthogonal to the tangent
struct SplineSection {
    float3 c;
    float3 t;
    float3 w;
};

// The B-spline's weights at u, and of its derivative: on control points
// p0..p3 (x = u + 1) before e, p1..p4 (x = u) from e on
struct SplineWeights {
    bool  lo;
    float b0, b1, b2, b3;
    float d0, d1, d2, d3;
    float e0, e1, e2, e3;       // of its second derivative
};

SplineWeights spline_weights(float u) {
    SplineWeights k;
    k.lo = u < 0.0;
    float x  = k.lo ? u + 1.0 : u;
    float x2 = x * x, x3 = x2 * x, ix = 1.0 - x;
    k.b0 = ix * ix * ix * (1.0 / 6.0);
    k.b1 = (3.0 * x3 - 6.0 * x2 + 4.0) * (1.0 / 6.0);
    k.b2 = (-3.0 * x3 + 3.0 * x2 + 3.0 * x + 1.0) * (1.0 / 6.0);
    k.b3 = x3 * (1.0 / 6.0);
    k.d0 = -0.5 * ix * ix;
    k.d1 = 1.5 * x2 - 2.0 * x;
    k.d2 = -1.5 * x2 + x + 0.5;
    k.d3 = 0.5 * x2;
    k.e0 = ix;
    k.e1 = 3.0 * x - 2.0;
    k.e2 = 1.0 - 3.0 * x;
    k.e3 = x;
    return k;
}

// The curve at u: the cross section's centre. spline_section's centre is
// this, so pieces tested by their centres alone meet the sections they get.
float3 spline_centre(SplinePart sp, float u) {
    SplineWeights k = spline_weights(u);
    float3 q0 = k.lo ? sp.p0 : sp.p1, q1 = k.lo ? sp.p1 : sp.p2, q2 = k.lo ? sp.p2 : sp.p3, q3 = k.lo ? sp.p3 : sp.p4;
    return q0 * k.b0 + q1 * k.b1 + q2 * k.b2 + q3 * k.b3;
}

SplineSection spline_section(SplinePart sp, float u) {
    SplineWeights k = spline_weights(u);
    float3 q0 = k.lo ? sp.p0 : sp.p1, q1 = k.lo ? sp.p1 : sp.p2, q2 = k.lo ? sp.p2 : sp.p3, q3 = k.lo ? sp.p3 : sp.p4;
    float3 v0 = k.lo ? sp.w0 : sp.w1, v1 = k.lo ? sp.w1 : sp.w2, v2 = k.lo ? sp.w2 : sp.w3, v3 = k.lo ? sp.w3 : sp.w4;
    SplineSection cs;
    cs.c = spline_centre(sp, u);
    cs.t = q0 * k.d0 + q1 * k.d1 + q2 * k.d2 + q3 * k.d3;
    float3 w  = v0 * k.b0 + v1 * k.b1 + v2 * k.b2 + v3 * k.b3;
    float  tt = dot(cs.t, cs.t);
    if (tt > 0.0) w -= cs.t * (dot(w, cs.t) / tt);
    float  ww = dot(w, w);
    if (!(ww > 1e-24)) {                                // no axis (tubes): any normal of the tangent
        float3 tn = tt > 0.0 ? cs.t / sqrt(tt) : float3(0.0, 0.0, 1.0);
        w = abs(tn.x) < 0.9 ? float3(1.0, 0.0, 0.0) : float3(0.0, 1.0, 0.0);
        w -= tn * dot(w, tn);
        ww = dot(w, w);
    }
    cs.w = w * rsqrt(ww);
    return cs;
}

// u at the end of piece k of n (k = 0: the part's start, n: its end); one
// expression wherever a piece's ends are needed
float spline_piece_u(SplinePart sp, uint k, uint n) {
    return lerp(sp.u0, sp.u1, float(k) / float(n));
}

// Pieces of a part (g: its bounds): as many as keep the straight pieces
// within about a fifth of a pixel of the curve. A piece of 1/K of the part
// strays about e / K^2 from it, e the part's own bend away from straight at
// its middle: the centre's, the rims' (the turn of the axis bends them
// further) and, on a ribbon, the faces' (the turn of the normal). Radii and
// widths are linear along the part, as along its pieces: no error. Straight
// parts are one piece, and exact.
uint spline_pieces(ViewData v, SplinePart sp, Geo g) {
    float px = M(v, proj).c1.y * 0.5 * v.viewport.y;
    if (v.ortho == 0u) px /= max(-0.5 * (g.a.z + g.b.z) - g.r, 1e-3);
    float  au = spline_width_rate(sp), am = 0.5 * (sp.a0 + sp.a1), rm = 0.5 * (sp.r0 + sp.r1);
    SplineFrame s0 = spline_frame(sp, sp.u0, sp.a0, au, sp.ru0);
    SplineFrame sm = spline_frame(sp, 0.5 * (sp.u0 + sp.u1), am, au, 0.5 * (sp.ru0 + sp.ru1));
    SplineFrame s1 = spline_frame(sp, sp.u1, sp.a1, au, sp.ru1);
    float3 dc = sm.c - 0.5 * (s0.c + s1.c);
    float3 dw = am * sm.w - 0.5 * (sp.a0 * s0.w + sp.a1 * s1.w);
    float  e  = max(length(dc + dw), length(dc - dw));
    if (sp.a0 > 0.0 || sp.a1 > 0.0) {                   // a tube's normal is arbitrary, and it has no faces
        float3 dfp = rm * normalize(sm.fp) - 0.5 * (sp.r0 * normalize(s0.fp) + sp.r1 * normalize(s1.fp));   // the faces' corners
        float3 dfm = rm * normalize(sm.fm) - 0.5 * (sp.r0 * normalize(s0.fm) + sp.r1 * normalize(s1.fm));
        float3 dbp = rm * normalize(sm.bp) - 0.5 * (sp.r0 * normalize(s0.bp) + sp.r1 * normalize(s1.bp));
        float3 dbm = rm * normalize(sm.bm) - 0.5 * (sp.r0 * normalize(s0.bm) + sp.r1 * normalize(s1.bm));
        e = max(e, max(max(length(dc + dw + dfp), length(dc + dw + dbp)), max(length(dc - dw + dfm), length(dc - dw + dbm))));
    }
    return uint(clamp(ceil(sqrt(e * px * 5.0)), 1.0, 16.0));
}

// The pieces of element e's part, element_geo's g: the rasters' count, as
// eval_hit counts them
uint spline_element_pieces(ViewData v, Batch b, uint inst, uint e, Geo g) {
    return spline_pieces(v, spline_part(v, b, inst, e, true), g);
}

// The ray against the bilinear patch q00 q10 q11 q01: u along q00 -> q10,
// v along q00 -> q01 (Reshetov, "Cool Patches", Ray Tracing Gems ch. 8).
// A piece's face is one: two triangles would stand off a twisting face by
// half its turn times its width, a patch by an eighth of the turn squared.
void patch_root(Ray r, float3 a0, float3 b0, float3 e00, float3 e11, float u, inout float t, inout float hu, inout float hv) {
    if (!(u >= 0.0 && u <= 1.0)) return;
    float3 pa = lerp(a0, b0, u);
    float3 pb = lerp(e00, e11, u);
    float3 n  = cross(r.d, pb);
    float  dd = dot(n, n);
    n = cross(n, pa);
    float tt = dot(n, pb) / dd;
    float vv = dot(n, r.d);
    if (tt > 0.0 && vv >= 0.0 && vv <= dd && (t < 0.0 || tt < t)) { t = tt; hu = u; hv = vv / dd; }
}

float patch_hit(Ray r, float3 q00, float3 q10, float3 q11, float3 q01, out float hu, out float hv) {
    hu = 0.0;
    hv = 0.0;
    float3 e10 = q10 - q00, e11 = q11 - q10, e00 = q01 - q00;
    float3 qn  = cross(e10, q01 - q11);
    float3 a0  = q00 - r.o, b0 = q10 - r.o;
    float  a   = dot(cross(a0, r.d), e00);              // a + b u + c u^2 = 0
    float  c   = dot(qn, r.d);
    float  b   = dot(cross(b0, r.d), e11) - (a + c);
    float  det = b * b - 4.0 * a * c;
    if (det < 0.0) return -1.0;
    det = sqrt(det);
    float u1, u2;
    if (c == 0.0) { u1 = -a / b; u2 = -1.0; }          // planar: one root (a NaN or inf misses)
    else          { u1 = 0.5 * (-b - (b >= 0.0 ? det : -det)); u2 = a / u1; u1 /= c; }
    float t = -1.0;
    patch_root(r, a0, b0, e00, e11, u1, t, hu, hv);
    patch_root(r, a0, b0, e00, e11, u2, t, hu, hv);
    return t;
}

// The normal of a capsule's surface at p
float3 capsule_normal(float3 p, float3 pa, float3 pb) {
    float3 n = p - lerp(pa, pb, capsule_param(p, pa, pb));
    float  l = length(n);
    return l > 0.0 ? n / l : float3(0.0, 0.0, 1.0);
}

struct SplineHit {
    float  t;                   // < 0: none
    float3 n;                   // view space, unit
};

void spline_take_capsule(Ray r, inout SplineHit h, float3 pa, float3 pb, float rad) {
    float t = capsule_hit(r, pa, pb, rad);
    if (t >= 0.0 && (h.t < 0.0 || t < h.t)) { h.t = t; h.n = capsule_normal(r.o + r.d * t, pa, pb); }
}

// A round cone: the balls of radius ra at pa to rb at pb, the radius linear
// between (their hull). Quilez's intersector, the body's root taken whatever
// the sign of its leading term (steep rays along a tapering cone); with one
// ball inside the other, the balls alone. Unit r.d.
struct ConeHit {
    float  t;                   // < 0: none
    float3 n;
};

ConeHit round_cone_hit(Ray r, float3 pa, float3 pb, float ra, float rb) {
    ConeHit h;
    h.t = -1.0;
    h.n = float3(0.0, 0.0, 1.0);
    float3 ba = pb - pa, oa = r.o - pa, ob = r.o - pb;
    float  rr = ra - rb;
    float  m0 = dot(ba, ba), m1 = dot(ba, oa), m2 = dot(ba, r.d), m3 = dot(r.d, oa);
    float  m5 = dot(oa, oa), m6 = dot(ob, r.d), m7 = dot(ob, ob);
    float  d2 = m0 - rr * rr;
    if (d2 > 0.0) {
        float k2 = d2 - m2 * m2;
        float k1 = d2 * m3 - m1 * m2 + m2 * rr * ra;
        float k0 = d2 * m5 - m1 * m1 + m1 * rr * ra * 2.0 - m0 * ra * ra;
        float hh = k1 * k1 - k0 * k2;
        if (hh >= 0.0 && abs(k2) > 1e-12 * m0) {
            float t = (-sqrt(hh) - k1) / k2;
            float y = m1 - ra * rr + t * m2;
            if (y > 0.0 && y < d2 && t >= 0.0) { h.t = t; h.n = normalize(d2 * (oa + t * r.d) - ba * y); return h; }
        }
    }
    float h1 = m3 * m3 - m5 + ra * ra;
    float h2 = m6 * m6 - m7 + rb * rb;
    if (h1 > 0.0) {
        float t = -m3 - sqrt(h1);
        if (t >= 0.0) { h.t = t; h.n = (oa + t * r.d) / ra; }
    }
    if (h2 > 0.0) {
        float t = -m6 - sqrt(h2);
        if (t >= 0.0 && (h.t < 0.0 || t < h.t)) { h.t = t; h.n = (ob + t * r.d) / rb; }
    }
    return h;
}

// A rim: the round cone along a piece's edge, between the rim's points pa
// and pb, whose tangents (per u) are ta and tb; du: the piece's length in
// u. The cone is the chord of the rim; its normal is taken from the rim
// itself, a cubic (Hermite) through the ends with their tangents: the
// direction from the rim's point at the hit's place along the piece,
// across its tangent there, tilted by the cone's slope. So shading follows
// the curve, not the chord: no facets at the joints, and no wobble along a
// thin rim, whose chord's sag would turn the normal by sag / radius (6 to
// 12 degrees on the test's tight bends, at a fifth of a pixel of sag).
//
// ends: 1 if pa, 2 if pb is a capped end of the part. There the rim's ball
// is the cap's and keeps its own normal (it and the cap are one sphere,
// which either may win by an ulp).
void spline_take_rim(Ray r, inout SplineHit h, float3 pa, float3 pb, float ra, float rb, float3 ta, float3 tb, float du, uint ends) {
    ConeHit c = round_cone_hit(r, pa, pb, ra, rb);
    if (c.t >= 0.0 && (h.t < 0.0 || c.t < h.t)) {
        float3 p  = r.o + r.d * c.t;
        float  k  = capsule_param(p, pa, pb);
        h.t = c.t;
        if ((k <= 0.0 && (ends & 1u) != 0u) || (k >= 1.0 && (ends & 2u) != 0u)) { h.n = c.n; return; }
        float  k2 = k * k, k3 = k2 * k;
        float3 ma = ta * du, mb = tb * du;
        float3 q  = pa * (2.0 * k3 - 3.0 * k2 + 1.0) + ma * (k3 - 2.0 * k2 + k) + pb * (3.0 * k2 - 2.0 * k3) + mb * (k3 - k2);
        float3 qd = (pa - pb) * (6.0 * k2 - 6.0 * k) + ma * (3.0 * k2 - 4.0 * k + 1.0) + mb * (3.0 * k2 - 2.0 * k);
        float  qq = dot(qd, qd);
        float3 tg = qq > 0.0 ? qd * rsqrt(qq) : float3(0.0);
        float3 m  = p - q;
        m -= tg * dot(m, tg);                           // across the rim
        float  mm = dot(m, m);
        float3 ba = pb - pa;
        float  bb = dot(ba, ba);
        float  sa = bb > 0.0 ? clamp((ra - rb) * rsqrt(bb), -1.0, 1.0) : 0.0;
        if (mm > 1e-30) h.n = m * (rsqrt(mm) * sqrt(1.0 - sa * sa)) + tg * sa;
        else            h.n = c.n;
    }
}

// A face: from the section at the piece's start (q00 on the +axis rim,
// q01 on the -axis one) to its end (q10, q11). Its normal is bilinear in
// the corners' unnormalised ones (SplineFrame mp, mm): the swept surface's
// is linear across the face before normalising, and normalised corners
// would weigh a bent face's edges alike (11 degrees off at the middle of
// an arrow's base)
void spline_take_face(Ray r, inout SplineHit h, float3 q00, float3 q01, float3 q10, float3 q11,
                      float3 n00, float3 n01, float3 n10, float3 n11) {
    float u, v;
    float t = patch_hit(r, q00, q10, q11, q01, u, v);
    if (t >= 0.0 && (h.t < 0.0 || t < h.t)) { h.t = t; h.n = normalize(lerp(lerp(n00, n01, v), lerp(n10, n11, v), u)); }
}

// The swept face's unit normal at an edge: m (across the rim's tangent tg
// and the axis w), tilted along the rim's direction by the radius's rate
float3 spline_edge_normal(float3 m, float3 tg, float3 w, float ru) {
    float3 n  = normalize(m);
    float3 e  = tg - w * dot(tg, w);                    // along the surface, across w
    float  el = length(e);
    if (ru == 0.0 || !(el > 0.0)) return n;
    float  mu = clamp(-ru / el, -0.99, 0.99);
    return n * sqrt(1.0 - mu * mu) + e * (mu / el);
}

// The full cross section, what a piece's rims and faces are built from: its
// centre and axis, and at its two edges (s = +1 and -1, the rims) the swept
// surface's normal and the rim's tangent.
//
// The face is the swept segment S(u, s) = c + s a w moved out by the radius
// along its own normal, cross(dS/du, w) = |t| n + s a cross(w', w), with
// w' = alpha t^ + beta n: the axis turning about the tangent (beta, the
// twist) tilts the normal across the face, by atan(a beta / |t|) at the
// edges (13 degrees on a helix ribbon). That normal is the rim tube's there
// (it is across the rim's tangent t +- (a' w + a w')): a face built on the
// edge normals meets the rims flush, in position and in shading. Built on
// the centre's normal it met them at a crease of that angle, which SSAO
// draws as a groove along every twisting rim.
struct SplineFrame {
    float3 c, w;
    float3 n;                   // the centre's normal, cross(t, w), unit
    float3 fp, fc, fm;          // the front face's normal at the + edge, the
    float3 bp, bc, bm;          // middle and the - edge, and the back face's,
                                // each of the length that weighs it in the
                                // face's bilinear blend
    float3 tp, tm;              // the + and - rim's tangent (not unit)
};

// a: the half width at u; au, ru: the half width's and the radius's rate
// along u. A radius that changes along the part tilts the normal along the
// surface by asin(-ru / speed): balls of a shrinking radius face forward.
// (The rims' cones carry the same tilt.)
SplineFrame spline_frame(SplinePart sp, float u, float a, float au, float ru) {
    SplineSection cs = spline_section(sp, u);
    SplineWeights k  = spline_weights(u);
    float3 q0 = k.lo ? sp.p0 : sp.p1, q1 = k.lo ? sp.p1 : sp.p2, q2 = k.lo ? sp.p2 : sp.p3, q3 = k.lo ? sp.p3 : sp.p4;
    float3 v0 = k.lo ? sp.w0 : sp.w1, v1 = k.lo ? sp.w1 : sp.w2, v2 = k.lo ? sp.w2 : sp.w3, v3 = k.lo ? sp.w3 : sp.w4;
    float3 t2 = q0 * k.e0 + q1 * k.e1 + q2 * k.e2 + q3 * k.e3;      // t'
    float3 vb = v0 * k.b0 + v1 * k.b1 + v2 * k.b2 + v3 * k.b3;      // the axes' blend
    float3 vd = v0 * k.d0 + v1 * k.d1 + v2 * k.d2 + v3 * k.d3;      // and its rate
    float  tt = dot(cs.t, cs.t);
    float  tl = sqrt(tt);
    float3 tn = tl > 0.0 ? cs.t / tl : float3(0.0, 0.0, 1.0);
    SplineFrame fr;
    fr.c = cs.c;
    fr.w = cs.w;
    fr.n = normalize(cross(tn, cs.w));
    // w = (vb - t (vb.t)/(t.t)) / |...|: alpha = w'.t^ = -(w.t')/|t| (w.t = 0
    // throughout), beta = w'.n = (vb'.n - (vb.t)/(t.t) t'.n) / |...|
    float  kv = tt > 0.0 ? dot(vb, cs.t) / tt : 0.0;
    float  wl = length(vb - cs.t * kv);
    float  alpha = 0.0, beta = 0.0;
    if (wl > 1e-12 && tl > 0.0) {                       // no axis (tubes): no turn
        alpha = -dot(cs.w, t2) / tl;
        beta  = (dot(vd, fr.n) - kv * dot(t2, fr.n)) / wl;
    }
    // where the bend folds the face (|t| - a alpha <= 0, a curve tighter
    // than the half width) its normal is held on the centre's side
    float3 mp = max(tl + a * alpha, 0.1 * tl) * fr.n - (a * beta) * tn;
    float3 mm = max(tl - a * alpha, 0.1 * tl) * fr.n + (a * beta) * tn;
    float3 wd = alpha * tn + beta * fr.n;               // w'
    fr.tp = cs.t + au * cs.w + a * wd;
    fr.tm = cs.t - au * cs.w - a * wd;
    // the tilt along u is the same on both sides: the back face's normal is
    // not the front's negated
    float lp = length(mp), lm = length(mm);
    fr.fp = spline_edge_normal( mp, fr.tp, cs.w, ru) * lp;
    fr.fm = spline_edge_normal( mm, fr.tm, cs.w, ru) * lm;
    fr.bp = spline_edge_normal(-mp, fr.tp, cs.w, ru) * lp;
    fr.bm = spline_edge_normal(-mm, fr.tm, cs.w, ru) * lm;
    fr.fc = spline_edge_normal( fr.n, cs.t, cs.w, ru) * tl;     // s = 0: |t| n
    fr.bc = spline_edge_normal(-fr.n, cs.t, cs.w, ru) * tl;
    return fr;
}

// The half width's rate along u: linear over the part
float spline_width_rate(SplinePart sp) {
    return sp.u1 > sp.u0 ? (sp.a1 - sp.a0) / (sp.u1 - sp.u0) : 0.0;
}

// The ray against the part's pieces. Two phases per piece: the centres of
// its ends (cheap) give its bounding sphere; only a piece the ray may hit
// gets its full sections, and the next piece reuses the one they share.
SplineHit spline_hit(Ray r, SplinePart sp, uint pieces) {
    SplineHit h;
    h.t = -1.0;
    h.n = float3(0.0, 0.0, 1.0);
    float3      c0 = spline_centre(sp, sp.u0);
    float       a0 = sp.a0, r0 = sp.r0;
    float       au = spline_width_rate(sp);
    float       du = (sp.u1 - sp.u0) / float(pieces);
    SplineFrame f0;
    f0.c = f0.w = f0.n = f0.fp = f0.fc = f0.fm = f0.bp = f0.bc = f0.bm = f0.tp = f0.tm = float3(0.0);
    bool        have0 = false;
    if ((sp.caps & 1u) != 0u) {
        f0 = spline_frame(sp, sp.u0, a0, au, sp.ru0); have0 = true;
        spline_take_capsule(r, h, f0.c - f0.w * a0, f0.c + f0.w * a0, r0);
    }
    for (uint k = 1u; k <= pieces; ++k) {
        float  u1 = spline_piece_u(sp, k, pieces);
        float3 c1 = spline_centre(sp, u1);
        float  a1 = lerp(sp.a0, sp.a1, float(k) / float(pieces));
        float  r1 = lerp(sp.r0, sp.r1, float(k) / float(pieces));
        // the piece's bounding sphere: its chord, widened by the profile
        float3 ctr = 0.5 * (c0 + c1);
        float  br  = 0.5 * length(c1 - c0) + max(a0, a1) + max(r0, r1);
        float3 oc  = r.o - ctr;
        float  bb  = dot(oc, r.d);
        float  hh  = bb * bb - (dot(oc, oc) - br * br);
        bool   have1 = false;
        SplineFrame f1;
        f1.c = f1.w = f1.n = f1.fp = f1.fc = f1.fm = f1.bp = f1.bc = f1.bm = f1.tp = f1.tm = float3(0.0);
        if (hh >= 0.0 && (h.t < 0.0 || -bb - sqrt(hh) < h.t)) {
            if (!have0) { f0 = spline_frame(sp, spline_piece_u(sp, k - 1u, pieces), a0, au, lerp(sp.ru0, sp.ru1, float(k - 1u) / float(pieces))); have0 = true; }
            f1 = spline_frame(sp, u1, a1, au, lerp(sp.ru0, sp.ru1, float(k) / float(pieces))); have1 = true;
            float3 p0 = f0.c + f0.w * a0, m0 = f0.c - f0.w * a0;      // the rims' ends
            float3 p1 = f1.c + f1.w * a1, m1 = f1.c - f1.w * a1;
            uint   ends = (k == 1u ? sp.caps & 1u : 0u) | (k == pieces ? sp.caps & 2u : 0u);
            spline_take_rim(r, h, p0, p1, r0, r1, f0.tp, f1.tp, du, ends);
            if (a0 > 0.0 || a1 > 0.0) {
                spline_take_rim(r, h, m0, m1, r0, r1, f0.tm, f1.tm, du, ends);
                // each face in two halves across, at its middle line: the
                // normal turns across the face (by the twist and the bend),
                // so its offset points lie on an arc of the radius, whose
                // chord in one patch stood in by r (1 - cos(turn / 2))
                float3 fq0 = f0.c + normalize(f0.fc) * r0, fq1 = f1.c + normalize(f1.fc) * r1;
                float3 bq0 = f0.c + normalize(f0.bc) * r0, bq1 = f1.c + normalize(f1.bc) * r1;
                float3 fp0 = p0 + normalize(f0.fp) * r0, fp1 = p1 + normalize(f1.fp) * r1;
                float3 fm0 = m0 + normalize(f0.fm) * r0, fm1 = m1 + normalize(f1.fm) * r1;
                float3 bp0 = p0 + normalize(f0.bp) * r0, bp1 = p1 + normalize(f1.bp) * r1;
                float3 bm0 = m0 + normalize(f0.bm) * r0, bm1 = m1 + normalize(f1.bm) * r1;
                spline_take_face(r, h, fp0, fq0, fp1, fq1, f0.fp, f0.fc, f1.fp, f1.fc);
                spline_take_face(r, h, fq0, fm0, fq1, fm1, f0.fc, f0.fm, f1.fc, f1.fm);
                spline_take_face(r, h, bp0, bq0, bp1, bq1, f0.bp, f0.bc, f1.bp, f1.bc);
                spline_take_face(r, h, bq0, bm0, bq1, bm1, f0.bc, f0.bm, f1.bc, f1.bm);
            }
        }
        c0 = c1; a0 = a1; r0 = r1; f0 = f1; have0 = have1;
    }
    if ((sp.caps & 2u) != 0u) {
        if (!have0) f0 = spline_frame(sp, spline_piece_u(sp, pieces, pieces), a0, au, sp.ru1);
        spline_take_capsule(r, h, f0.c - f0.w * a0, f0.c + f0.w * a0, r0);
    }
    return h;
}

float seg_dist(float3 p, float3 a, float3 ab, float ll) {
    float3 pq = p - a;
    return length(pq - ab * saturate(dot(pq, ab) / ll));
}

// The part's curve lies in the hull of each of its halves' Bezier control
// points: the ends (the curve at u0, 0 and u1) and two inner points each.
// Far tighter than the B-spline's own control points, which on a helix
// swing out two turns' worth (4.2 A off the chord where the curve is 0.5).
struct SplineHull {
    float3 ca, cb;              // the curve at u0 and u1: the chord
    float3 m;                   // the curve at u = 0 (the control point's own)
    float3 l1, l2;              // before e, u in [-0.5, 0]
    float3 h1, h2;              // from e on, u in [0, 0.5]
    bool   lo, hi;              // which halves the part has
};

SplineHull spline_hull(float3 p0, float3 p1, float3 p2, float3 p3, float3 p4, bool lo, bool hi) {
    SplineHull k;
    k.lo = lo; k.hi = hi;
    k.m  = (p1 + 4.0 * p2 + p3) * (1.0 / 6.0);
    k.l1 = (9.0 * p1 + 14.0 * p2 + p3) * (1.0 / 24.0);
    k.l2 = (3.0 * p1 + 8.0 * p2 + p3) * (1.0 / 12.0);
    k.h1 = (p1 + 8.0 * p2 + 3.0 * p3) * (1.0 / 12.0);
    k.h2 = (p1 + 14.0 * p2 + 9.0 * p3) * (1.0 / 24.0);
    k.ca = lo ? (p0 + 23.0 * p1 + 23.0 * p2 + p3) * (1.0 / 48.0) : k.m;
    k.cb = hi ? (p1 + 23.0 * p2 + 23.0 * p3 + p4) * (1.0 / 48.0) : k.m;
    return k;
}

// A bounding capsule of the part, in view space: its curve's ends (as the
// pieces have them), widened by the farthest of its hull and the profile
void spline_bounds(SplinePart sp, out float3 a, out float3 bnd_b, out float rad) {
    a = spline_centre(sp, sp.u0);
    bnd_b = spline_centre(sp, sp.u1);
    SplineHull k = spline_hull(sp.p0, sp.p1, sp.p2, sp.p3, sp.p4, sp.u0 < 0.0, sp.u1 > 0.0);
    float3 ab = bnd_b - a;
    float  ll = max(dot(ab, ab), 1e-30);
    float  d = seg_dist(k.m, a, ab, ll);
    if (k.lo) d = max(d, max(seg_dist(k.l1, a, ab, ll), seg_dist(k.l2, a, ab, ll)));
    if (k.hi) d = max(d, max(seg_dist(k.h1, a, ab, ll), seg_dist(k.h2, a, ab, ll)));
    rad = d + max(sp.a0, sp.a1) + max(sp.r0, sp.r1);
}

// Spheres and capsules (spline parts: spline_geo). With `fade` false (the
// state lists, which do not depend on the camera) the batch's fade is not
// applied.
Geo element_geo(ViewData v, Batch b, uint inst, uint e, bool fade) {
    Geo g;
    g.ok   = false;
    g.faded = false;
    g.kind = b.kind == KIND_SPHERES ? 0u : 1u;
    g.st   = 0u;
    g.a    = float3(0.0);
    g.b    = float3(0.0);
    g.r    = 0.0;
    Mat   T = load_transform(b.transforms, inst);
    float s = length(T.c0.xyz);
    if (b.kind == KIND_SPHERES) {
        uint pi = point_of(b, e, 0u, 1u);
        uint ai = attr_of(b, pi);
        if (!attr_visible(b, ai)) return g;
        g.st = attr_state(b, ai);
        g.r  = sphere_radius(b, pi) * s;
        g.a  = mulp(M(v, view), mulp(T, load_pos(b, pi)));
        g.b  = g.a;
        if (fade && (b.flags & F_FADE) != 0u && g.r > 0.0) {
            float2 t = fade_ts(b, -g.a.z);
            if ((b.flags & F_FADE_SHRINK) != 0u) g.r *= fade_size(t);
            g.faded = g.r <= 0.0 || ((b.flags & F_FADE_THIN) != 0u && !fade_keeps(t, ai));
        }
    } else if (b.kind == KIND_CAPSULES) {
        uint p0 = point_of(b, e, 0u, 2u);
        uint p1 = point_of(b, e, 1u, 2u);
        uint a0 = attr_of(b, p0);
        uint a1 = attr_of(b, p1);
        if (!attr_visible(b, a0) || !attr_visible(b, a1)) return g;
        g.st = attr_state(b, a0) & attr_state(b, a1);
        g.r  = capsule_radius(b, e) * s;
        g.a  = mulp(M(v, view), mulp(T, load_pos(b, p0)));
        g.b  = mulp(M(v, view), mulp(T, load_pos(b, p1)));
        if (fade && (b.flags & F_FADE) != 0u && g.r > 0.0) {
            float2 t = fade_ts(b, -0.5 * (g.a.z + g.b.z));      // the midpoint
            if ((b.flags & F_FADE_SHRINK) != 0u) g.r *= fade_size(t);
            g.faded = g.r <= 0.0 || ((b.flags & F_FADE_THIN) != 0u && !(fade_keeps(t, a0) && fade_keeps(t, a1)));
        }
    } else {
        return g;
    }
    g.ok = g.r > 0.0 && !g.faded;
    return g;
}

// A spline part's bounds: a capsule around it (the box proxy's), as
// element_geo's for capsules
Geo spline_geo(ViewData v, Batch b, uint inst, uint e, bool fade) {
    Geo g;
    g.ok    = false;
    g.faded = false;
    g.kind  = 2u;
    g.st    = 0u;
    g.a     = float3(0.0);
    g.b     = float3(0.0);
    g.r     = 0.0;
    uint ie = b.index[e];
    uint ai = attr_of(b, ie & SPLINE_POINT);
    if (!attr_visible(b, ai)) return g;
    g.st = attr_state(b, ai);
    if ((ie & (SPLINE_FIRST | SPLINE_LAST)) == (SPLINE_FIRST | SPLINE_LAST)) return g;   // alone
    SplinePart sp = spline_part(v, b, inst, e, fade);
    if (fade && (b.flags & F_FADE) != 0u) {
        float2 t = fade_ts(b, -sp.p2.z);
        g.faded = ((b.flags & F_FADE_SHRINK) != 0u && fade_size(t) <= 0.0) || ((b.flags & F_FADE_THIN) != 0u && !fade_keeps(t, ai));
    }
    if (!sp.ok || g.faded) return g;
    spline_bounds(sp, g.a, g.b, g.r);
    g.ok = g.r > 0.0 && !g.faded;
    return g;
}

// ---------------------------------------------------------------------------
// Polygons (md_render_polygon_style_t): a flat convex polygon of n corners
// (3 .. 8, the same for a whole batch), made solid by a ball of the radius.
// Its plane: Newell's normal of the corners about their centroid (the fit of
// a polygon given in order, oriented by its winding); the corners projected
// onto it. The solid is the plate, the polygon moved out by the radius to
// either side, and a capsule along each edge: the Minkowski sum, for a
// convex polygon. Every stage builds the polygon with poly_part, so they
// agree to the bit. The corners live in registers: every index into q is a
// constant of a loop Slang unrolls ([unroll] is only a hint to the driver:
// lavapipe's pick kernel crashed on the array indexed at run time, and
// elsewhere it would have been local memory).
// ---------------------------------------------------------------------------

static const uint POLY_MAX = 8u;

struct PolyPart {
    bool   ok;
    bool   plate;               // a plane: corners not on a line
    uint   n;
    float3 q[POLY_MAX];         // projected onto the plane, scaled about c (shrink)
    float3 c;                   // centroid
    float3 nrm;                 // unit normal, by the winding
    float  r;
};

uint poly_n(Batch b) { return clamp(b.poly_n, 3u, POLY_MAX); }

// The corner after k (constant k)
float3 poly_next(PolyPart pp, uint k) {
    return k + 1u < pp.n ? pp.q[min(k + 1u, POLY_MAX - 1u)] : pp.q[0];
}

float3 poly_any_perpendicular(float3 t) {
    float3 v = abs(t.x) < 0.9 ? float3(1.0, 0.0, 0.0) : float3(0.0, 1.0, 0.0);
    return normalize(cross(t, v));
}

// Centroid, plane, projection: of corners q[0 .. n) in any space (a
// similarity maps the plane of one space to the plane of the other)
void poly_plane(inout PolyPart pp) {
    float3 c = float3(0.0);
    [ForceUnroll] for (uint k = 0u; k < POLY_MAX; ++k) if (k < pp.n) c += pp.q[k];
    c /= float(pp.n);
    float3 nv = float3(0.0);
    float  sz = 0.0;
    [ForceUnroll] for (uint k = 0u; k < POLY_MAX; ++k) {
        if (k < pp.n) {
            nv += cross(pp.q[k] - c, poly_next(pp, k) - c);
            sz = max(sz, length(pp.q[k] - c));
        }
    }
    float nl = length(nv);
    pp.plate = nl > 1e-6 * sz * sz;
    if (pp.plate) {
        pp.nrm = nv / nl;
        [ForceUnroll] for (uint k = 0u; k < POLY_MAX; ++k) if (k < pp.n) pp.q[k] -= pp.nrm * dot(pp.q[k] - c, pp.nrm);
    } else {
        // on a line (or a point): the capsules of its edges alone
        float3 d = pp.q[0] - c;
        pp.nrm = poly_any_perpendicular(dot(d, d) > 0.0 ? normalize(d) : float3(0.0, 0.0, 1.0));
    }
    pp.c = c;
}

PolyPart poly_part(ViewData v, Batch b, uint inst, uint e, bool fade) {
    PolyPart pp;
    pp.ok = false; pp.plate = false; pp.r = 0.0;
    pp.n  = poly_n(b);
    Mat T = load_transform(b.transforms, inst);
    Mat V = M(v, view);
    [ForceUnroll] for (uint k = 0u; k < POLY_MAX; ++k) {
        pp.q[k] = float3(0.0);
        if (k < pp.n) pp.q[k] = mulp(V, mulp(T, load_pos(b, point_of(b, e, k, pp.n))));
    }
    poly_plane(pp);
    float f = 1.0;
    if (fade && (b.flags & F_FADE) != 0u && (b.flags & F_FADE_SHRINK) != 0u) f = fade_size(fade_ts(b, -pp.c.z));
    if (f != 1.0) {
        [ForceUnroll] for (uint k = 0u; k < POLY_MAX; ++k) if (k < pp.n) pp.q[k] = pp.c + (pp.q[k] - pp.c) * f;
    }
    pp.r  = b.s0 * length(T.c0.xyz) * f;
    pp.ok = pp.r > 0.0;
    return pp;
}

// Visible when all its corners are, its state bits all theirs
bool poly_visible(Batch b, uint e, out uint st) {
    uint n = poly_n(b);
    bool vis = true;
    st = 0xFFFFFFFFu;
    for (uint k = 0u; k < n; ++k) {
        uint ai = attr_of(b, point_of(b, e, k, n));
        vis = vis && attr_visible(b, ai);
        st &= attr_state(b, ai);
    }
    return vis;
}

// The ray against the polygon: the plate's faces, then a capsule per edge
SplineHit poly_hit(Ray r, PolyPart pp) {
    SplineHit h;
    h.t = -1.0;
    h.n = float3(0.0, 0.0, 1.0);
    float dn = dot(r.d, pp.nrm);
    if (pp.plate && abs(dn) > 1e-12) {
        float sg = dn < 0.0 ? 1.0 : -1.0;           // the face towards the ray
        float t  = dot(pp.c + pp.nrm * (pp.r * sg) - r.o, pp.nrm) / dn;
        if (t >= 0.0) {
            float3 p = r.o + r.d * t - pp.nrm * (pp.r * sg);     // onto the polygon's plane
            bool inside = true;
            [ForceUnroll] for (uint k = 0u; k < POLY_MAX; ++k) {
                if (k < pp.n) inside = inside && dot(cross(poly_next(pp, k) - pp.q[k], p - pp.q[k]), pp.nrm) >= 0.0;
            }
            if (inside) { h.t = t; h.n = pp.nrm * sg; }
        }
    }
    if (h.t >= 0.0) return h;                       // a face, where it is, is nearer than any rim
    [ForceUnroll] for (uint k = 0u; k < POLY_MAX; ++k) {
        if (k < pp.n) spline_take_capsule(r, h, pp.q[k], poly_next(pp, k), pp.r);
    }
    return h;
}

// Its bounding sphere: as a capsule of length 0
Geo poly_geo_part(ViewData v, Batch b, uint inst, uint e, bool fade, out PolyPart pp) {
    Geo g;
    g.ok = false; g.faded = false; g.kind = 3u; g.st = 0u; g.a = float3(0.0); g.b = float3(0.0); g.r = 0.0;
    pp = poly_part(v, b, inst, e, fade);
    uint st;
    if (!poly_visible(b, e, st)) return g;
    g.st = st;
    if (fade && (b.flags & F_FADE) != 0u) {
        float2 t = fade_ts(b, -pp.c.z);
        bool keep = true;
        if ((b.flags & F_FADE_THIN) != 0u) {
            for (uint k = 0u; k < pp.n; ++k) keep = keep && fade_keeps(t, attr_of(b, point_of(b, e, k, pp.n)));
        }
        g.faded = ((b.flags & F_FADE_SHRINK) != 0u && fade_size(t) <= 0.0) || !keep;
    }
    if (!pp.ok || g.faded) return g;
    float d = 0.0;
    [ForceUnroll] for (uint k = 0u; k < POLY_MAX; ++k) if (k < pp.n) d = max(d, length(pp.q[k] - pp.c));
    g.a  = pp.c;
    g.b  = pp.c;
    g.r  = d + pp.r;
    g.ok = true;
    return g;
}

Geo poly_geo(ViewData v, Batch b, uint inst, uint e, bool fade) {
    PolyPart pp;
    return poly_geo_part(v, b, inst, e, fade, pp);
}

// Corner c of the box around the polygon, in clip space: in its plane the
// rectangle of its corners along u (towards the first) and w = n x u,
// widened by the radius; across it, the radius. (u, w, n) is right-handed,
// as the capsules' box frame: the same indices face out.
float4 poly_box_corner(ViewData v, PolyPart pp, uint c) {
    float3 u  = pp.q[0] - pp.c;
    u -= pp.nrm * dot(u, pp.nrm);
    float  ul = length(u);
    u = ul > 1e-20 ? u / ul : poly_any_perpendicular(pp.nrm);
    float3 w  = cross(pp.nrm, u);
    float2 lo = float2(0.0), hi = float2(0.0);     // the centroid is inside
    [ForceUnroll] for (uint k = 0u; k < POLY_MAX; ++k) {
        if (k < pp.n) {
            float2 x = float2(dot(pp.q[k] - pp.c, u), dot(pp.q[k] - pp.c, w));
            lo = min(lo, x); hi = max(hi, x);
        }
    }
    float  br = pp.r * 1.002;               // see sphere_vs: not just touching
    float3 p  = pp.c
              + u * (((c & 1u) != 0u) ? hi.x + br : lo.x - br)
              + w * (((c & 2u) != 0u) ? hi.y + br : lo.y - br)
              + pp.nrm * (((c & 4u) != 0u) ? br : -br);
    return mul4(M(v, proj), float4(p, 1.0));
}

// Screen bounds of a view-space sphere (the bounding sphere of an element),
// from the unjittered projection, padded by a pixel for the jitter.
struct Bounds {
    bool   outside;             // wholly outside the frustum
    bool   near_clip;           // crosses the near plane: no bounds, no occlusion test
    float2 pmin;                // pixels, top-left origin
    float2 pmax;
    float  dnear;               // depth of its nearest point
    float  dnear_m;             // ... of a point occl_margin radii nearer still (phase 1's test)
};

Bounds sphere_bounds(ViewData v, float3 c, float r) {
    Bounds bd;
    bd.outside   = false;
    bd.near_clip = false;
    bd.pmin      = float2(0.0);
    bd.pmax      = v.viewport.xy;
    bd.dnear     = 1.0;
    bd.dnear_m   = 1.0;
    if (c.z - r > v.near_z) { bd.outside = true; return bd; }
    if (c.z + r > v.near_z) { bd.near_clip = true; return bd; }
    Mat P = M(v, proj_cur);
    float2 lo, hi;
    if (v.ortho != 0u) {
        float2 nc = float2(P.c0.x * c.x, P.c1.y * c.y) + P.c3.xy;
        float2 nr = abs(float2(P.c0.x, P.c1.y)) * r;
        lo = nc - nr;
        hi = nc + nr;
    } else {
        // The tangents from the eye in the xz and yz planes: x / -z of the two
        // tangent points, exact for a sphere in front of the near plane.
        float zf = -c.z;
        float2 cc = c.xy;
        float2 d2 = cc * cc + zf * zf;
        float2 t  = sqrt(max(d2 - r * r, 0.0));
        float2 s0 = (t * cc - r * zf) / (t * zf + r * cc);
        float2 s1 = (t * cc + r * zf) / (t * zf - r * cc);
        float2 smin = min(s0, s1);
        float2 smax = max(s0, s1);
        float2 sc = float2(P.c0.x, P.c1.y);
        float2 off = -P.c2.xy;              // off-centre projections
        float2 a0 = smin * sc + off;
        float2 a1 = smax * sc + off;
        lo = min(a0, a1);
        hi = max(a0, a1);
    }
    // ndc -> pixels (y flips)
    float2 wh = v.viewport.xy;
    bd.pmin = float2((lo.x * 0.5 + 0.5) * wh.x, (0.5 - hi.y * 0.5) * wh.y) - 1.0;
    bd.pmax = float2((hi.x * 0.5 + 0.5) * wh.x, (0.5 - lo.y * 0.5) * wh.y) + 1.0;
    float4 cz = mul4(P, float4(c.x, c.y, c.z + r, 1.0));
    bd.dnear = cz.z / cz.w;
    bd.dnear_m = bd.dnear;
    if (v.occl_margin > 0.0) {
        float zm = c.z + r * (1.0 + v.occl_margin);
        float4 cm = mul4(P, float4(c.x, c.y, zm, 1.0));
        bd.dnear_m = zm < v.near_z ? cm.z / cm.w : 2.0;     // past the near plane: never occluded
    }
    if (bd.pmax.x <= 0.0 || bd.pmax.y <= 0.0 || bd.pmin.x >= wh.x || bd.pmin.y >= wh.y || bd.dnear < 0.0) bd.outside = true;
    return bd;
}

// A capsule is the convex hull of its two end spheres, so its rectangle is
// the union of theirs (not the bounding sphere's: for a stick that is up to
// (length + 2r)^2 against (length + 2r) x 2r)
Bounds geo_bounds(ViewData v, Geo g) {
    if (g.kind == 0u) return sphere_bounds(v, g.a, g.r);
    Bounds bd;
    bd.outside   = false;
    bd.near_clip = false;
    bd.pmin      = float2(0.0);
    bd.pmax      = v.viewport.xy;
    bd.dnear     = 1.0;
    bd.dnear_m   = 1.0;
    if (min(g.a.z, g.b.z) - g.r > v.near_z) { bd.outside = true; return bd; }
    if (max(g.a.z, g.b.z) + g.r > v.near_z) { bd.near_clip = true; return bd; }
    Bounds b0 = sphere_bounds(v, g.a, g.r);
    Bounds b1 = sphere_bounds(v, g.b, g.r);
    bd.pmin  = min(b0.pmin, b1.pmin);
    bd.pmax  = max(b0.pmax, b1.pmax);
    bd.dnear   = max(b0.dnear, b1.dnear);
    bd.dnear_m = max(b0.dnear_m, b1.dnear_m);
    float2 wh = v.viewport.xy;
    if (bd.pmax.x <= 0.0 || bd.pmax.y <= 0.0 || bd.pmin.x >= wh.x || bd.pmin.y >= wh.y || bd.dnear < 0.0) bd.outside = true;
    return bd;
}

// Occluded by the previous raster: nearer than nothing in its rectangle. The
// pyramid holds the farthest (min, under reverse Z) depth of each texel's
// pixels; level L texels are 2^(L+1) pixels across, and the level is the one
// at which the rectangle spans at most 2 x 2 texels.
bool hiz_occluded(ViewData v, DescriptorHandle<Texture2D<float>> hiz, Bounds bd, bool margin) {
    if (v.hiz_valid == 0u || bd.near_clip) return false;
    float2 lo = max(bd.pmin, float2(0.0));
    float2 hi = min(bd.pmax, v.viewport.xy);
    float  w  = max(hi.x - lo.x, hi.y - lo.y);
    uint   L  = w <= 2.0 ? 0u : uint(ceil(log2(w))) - 1u;
    L = min(L, v.hiz_levels - 1u);
    uint dw = v.hiz_w, dh = v.hiz_h;
    for (uint j = 0u; j < L; ++j) { dw = max(1u, (dw + 1u) >> 1); dh = max(1u, (dh + 1u) >> 1); }
    float s  = float(1u << (L + 1u));
    int2  t0 = min(int2(floor(lo / s)), int2(int(dw) - 1, int(dh) - 1));
    int2  t1 = min(int2(floor((hi - 1e-3) / s)), int2(int(dw) - 1, int(dh) - 1));
    t1 = max(t0, t1);
    float m = hiz.Load(int3(t0.x, t0.y, int(L)));
    m = min(m, hiz.Load(int3(t1.x, t0.y, int(L))));
    m = min(m, hiz.Load(int3(t0.x, t1.y, int(L))));
    m = min(m, hiz.Load(int3(t1.x, t1.y, int(L))));
    return (margin ? bd.dnear_m : bd.dnear) < m;
}

// The pixel centres an element can reach: its exact rectangle (the bounds
// without their pad) widened by the jitter this raster applies and a hair
// for rounding. x0, y0, x1, y1, inclusive; none when x1 < x0 or y1 < y0.
int4 reach_px(ViewData v, Bounds bd) {
    float  g  = v.jitter_px + 0.05;
    float2 lo = bd.pmin + (1.0 - g);
    float2 hi = bd.pmax - (1.0 - g);
    return int4(int2(ceil(lo - 0.5)), int2(floor(hi - 0.5)));
}

// ---------------------------------------------------------------------------
// Cull
//
// Two levels. Clusters are tested first (frustum, depth pyramid); the
// elements of the clusters that pass are tested one by one and sorted into
// the hardware, small-primitive and occluded lists. Each phase of occlusion
// does both levels:
//
//   phase 1   clusters against the previous raster's pyramid; elements of
//             the visible ones likewise
//   phase 2   clusters occluded in phase 1, and elements occluded in phase 1,
//             against this raster's pyramid
//   state     clusters with a state bit, then their elements with one, into
//             the state lists for the masks: occluded or not, and not tested
//             against the frustum (the mask raster clips), so the lists hold
//             for any camera and are rebuilt only when the state bits or the
//             geometry change
//
// Statistics count elements, each exactly once: elements of hidden clusters
// and hidden elements of visible ones at the cluster level, the rest where
// their fate is decided.
// ---------------------------------------------------------------------------

static const uint OUT_HIDDEN   = 0u;
static const uint OUT_OUTSIDE  = 1u;
static const uint OUT_OCCLUDED = 2u;
static const uint OUT_HW       = 3u;
static const uint OUT_SMALL    = 4u;
static const uint OUT_NONE     = 5u;   // a lane without an element
static const uint OUT_LOD      = 6u;   // dropped by its batch's fade


// Appends with one atomic per wave: the index for this lane (valid when it wants one)
uint wave_append(uint* counter, bool want) {
    uint n    = WaveActiveCountBits(want);
    uint base = 0u;
    if (WaveIsFirstLane() && n > 0u) base = ((Atomic<uint>*)counter)->add(n);
    base = WaveReadLaneFirst(base);
    return base + WavePrefixSum(want ? 1u : 0u);    // WavePrefixCountBits has no Metal lowering
}

void wave_add(uint* counter, uint value) {
    uint n = WaveActiveSum(value);
    if (WaveIsFirstLane() && n > 0u) ((Atomic<uint>*)counter)->add(n);
}

void wave_count(uint* counter, bool want) { wave_add(counter, want ? 1u : 0u); }

uint grid_threads(uint n) { return 64u * min((n + 63u) / 64u, 65535u); }
uint grid_groups(uint n)  { return min(n, 65535u); }
uint cluster_groups(uint n) { return min((n + 63u) / 64u, 65535u); }   // the element cull: 64 clusters a group

// The batch of work item i: the last with offset <= i (batches without work
// share their offset with the next one, which wins)
uint find_batch_work(Batch* batches, uint count, uint i) {
    uint lo = 0u, hi = count;
    while (hi - lo > 1u) { uint mid = (lo + hi) >> 1; if (batches[mid].work_offset <= i) lo = mid; else hi = mid; }
    return lo;
}
// An element list entry (a global work index) as batch, instance, element
Elem decode_work(Batch* batches, uint batch_count, uint w) {
    Elem el;
    el.slot = find_batch_work(batches, batch_count, w);
    uint local = w - batches[el.slot].work_offset;
    uint cnt   = max(batches[el.slot].count, 1u);
    el.inst = local / cnt;
    el.e    = local - el.inst * cnt;
    return el;
}

// The visibility id of an element
uint2 elem_id(Batch* batches, Elem el) {
    return encode_id(el.slot, 0u, uint64_t(el.inst) * uint64_t(batches[el.slot].count) + uint64_t(el.e));
}

uint find_batch_cluster_work(Batch* batches, uint count, uint i) {
    uint lo = 0u, hi = count;
    while (hi - lo > 1u) { uint mid = (lo + hi) >> 1; if (batches[mid].cluster_work_offset <= i) lo = mid; else hi = mid; }
    return lo;
}
uint find_batch_cluster(Batch* batches, uint count, uint i) {
    uint lo = 0u, hi = count;
    while (hi - lo > 1u) { uint mid = (lo + hi) >> 1; if (batches[mid].cluster_base <= i) lo = mid; else hi = mid; }
    return lo;
}

uint2 cluster_range(Batch b, uint c) {
    uint cg = gated(b, F_CLUSTERS, c);
    if ((b.flags & F_CLUSTERS) != 0u) return uint2(b.cluster_offset[cg], b.cluster_offset[cg + 1u]);
    uint beg = c * 64u;
    return uint2(beg, min(beg + 64u, b.count));
}

// An element in the batch's own space: the visibility test, the state bits and
// the radius factor its style multiplies (spheres: the point radius, capsules:
// the per-capsule scale)
struct LocalElem {
    bool   ok;
    uint   st;
    float3 a;
    float3 b;
    float  factor;
};

// Without `geo` (a refit of the state bits alone) the positions are not read:
// over all elements, that is most of the traffic.
LocalElem element_local(Batch b, uint e, bool geo) {
    LocalElem l;
    l.ok = false; l.st = 0u; l.a = float3(0.0); l.b = float3(0.0); l.factor = 0.0;
    if (b.kind == KIND_SPHERES) {
        uint pi = point_of(b, e, 0u, 1u);
        uint ai = attr_of(b, pi);
        if (!attr_visible(b, ai)) return l;
        l.st = attr_state(b, ai);
        l.factor = point_radius(b, pi);
        if (geo) { l.a = load_pos(b, pi); l.b = l.a; }
        // From the factor, not from a second point_radius: unless the
        // compiler merges the two, that is a second load and, with a table, a
        // second lookup. The GTX 1060's refit at 10M atoms, a table against a
        // float radius: 0.24 ms with the second call (1c494694), 0.04 without
        // (66782714, a build with other differences too)
        l.ok = l.factor * b.s0 + b.s1 > 0.0;
    } else if (b.kind == KIND_CAPSULES) {
        uint p0 = point_of(b, e, 0u, 2u);
        uint p1 = point_of(b, e, 1u, 2u);
        uint a0 = attr_of(b, p0);
        uint a1 = attr_of(b, p1);
        if (!attr_visible(b, a0) || !attr_visible(b, a1)) return l;
        l.st = attr_state(b, a0) & attr_state(b, a1);
        l.factor = ((b.flags & F_RADIUS_SCALE) != 0u) ? b.radius_scale[gated(b, F_RADIUS_SCALE, e)] : 1.0;
        if (geo) { l.a = load_pos(b, p0); l.b = load_pos(b, p1); }
        l.ok = capsule_radius(b, e) > 0.0;
    }
    return l;
}

// A spline part's: the box of its curve's hull
LocalElem spline_local(Batch b, uint e, bool geo) {
    LocalElem l;
    l.ok = false; l.st = 0u; l.a = float3(0.0); l.b = float3(0.0); l.factor = 0.0;
    // a and b: the box of the part's curve, batch space; the factor its
    // widest half width and radius (the cluster's pad)
    uint ie = b.index[e];
    uint ai = attr_of(b, ie & SPLINE_POINT);
    if (!attr_visible(b, ai)) return l;
    bool hp = (ie & SPLINE_FIRST) == 0u, hn = (ie & SPLINE_LAST) == 0u;
    if (!hp && !hn) return l;
    l.st = attr_state(b, ai);
    float4 pf = spline_profile_of(b, e);
    float  pr = max(max(pf.z, pf.w), 0.0);
    l.factor = max(max(pf.x, pf.y), 0.0) + pr;
    if (geo) {
        SplineCtl c = spline_controls(b, e, ie);
        // the box of its hull (spline_hull): what the pieces stay in
        SplineHull k = spline_hull(c.p0, c.p1, c.p2, c.p3, c.p4, hp, hn);
        l.a = min(min(k.ca, k.cb), k.m);
        l.b = max(max(k.ca, k.cb), k.m);
        if (hp) { l.a = min(l.a, min(k.l1, k.l2)); l.b = max(l.b, max(k.l1, k.l2)); }
        if (hn) { l.a = min(l.a, min(k.h1, k.h2)); l.b = max(l.b, max(k.h1, k.h2)); }
    }
    l.ok = pr > 0.0;
    return l;
}

// A polygon's: the box of its corners on its plane (factor 1: the pad is
// the radius)
LocalElem poly_local(Batch b, uint e, bool geo) {
    LocalElem l;
    l.ok = false; l.st = 0u; l.a = float3(0.0); l.b = float3(0.0); l.factor = 1.0;
    uint st;
    if (!poly_visible(b, e, st)) return l;
    l.st = st;
    if (geo) {
        PolyPart pp;
        pp.ok = false; pp.plate = false; pp.r = 0.0;
        pp.n  = poly_n(b);
        [ForceUnroll] for (uint k = 0u; k < POLY_MAX; ++k) {
            pp.q[k] = float3(0.0);
            if (k < pp.n) pp.q[k] = load_pos(b, point_of(b, e, k, pp.n));
        }
        poly_plane(pp);
        l.a = pp.q[0];
        l.b = pp.q[0];
        [ForceUnroll] for (uint k = 1u; k < POLY_MAX; ++k) if (k < pp.n) { l.a = min(l.a, pp.q[k]); l.b = max(l.b, pp.q[k]); }
    }
    l.ok = b.s0 > 0.0;
    return l;
}

// Bounds of each cluster's visible elements, and its state bits. Eight
// threads per cluster, eight clusters per group: a thread takes every eighth
// element (several independent loads in flight), eight adjacent elements
// are one line, and the reduction's barriers serve eight clusters at once.
// (A thread per cluster made each load of a wave touch as many clusters, 1
// KB apart: on Pascal several times the bandwidth's time. A group per
// cluster left a load per thread between six barriers: latency.) Min, max,
// or and sum: the result does not depend on the order. Runs when the
// geometry changed (refit_mode 0), or when only the state bits did
// (refit_mode 1). Batches sharing bounds with an earlier one are skipped.
// Dispatched as refit_groups(clusters).
static const uint RF_LANES = 8u;
uint refit_groups(uint n) { return min((n + 7u) / 8u, 65535u); }

groupshared float4 rf_lo[64];       // xyz: min, w: the largest radius factor
groupshared float4 rf_hi[64];       // xyz: max
groupshared uint2  rf_cs[64];       // state bits, visible elements

void refit_body<let SPL : int>(uint3 gid, uint gi) {
    Args     a = MD_ARGS;
    ViewData v = *a.view;
    uint n    = v.total_clusters;
    bool geo  = v.refit_mode == 0u;
    uint lane = gi & (RF_LANES - 1u);
    for (uint base = gid.x * 8u; base < n; base += refit_groups(n) * 8u) {   // uniform over the group
        uint   i     = base + gi / RF_LANES;
        bool   valid = i < n;
        uint   ic    = valid ? i : n - 1u;                 // a cluster to read either way
        uint   slot  = find_batch_cluster(a.batches, v.batch_count, ic);
        Batch  b     = a.batches[slot];
        uint   c     = ic - b.cluster_base;
        valid = valid && (b.flags & F_SHARED) == 0u && kind_taken<SPL>(b.kind);
        float3 lo = float3(3.0e38), hi = float3(-3.0e38);
        float  fmax = 0.0;
        uint   count = 0u, st = 0u;
        if (valid) {
            uint2 rng = cluster_range(b, c);
            for (uint e = rng.x + lane; e < rng.y; e += RF_LANES) {
                LocalElem l;
                if (SPL != 0) l = b.kind == KIND_POLYGONS ? poly_local(b, e, geo) : spline_local(b, e, geo);
                else          l = element_local(b, e, geo);
                if (!l.ok) continue;
                st |= l.st;
                if (geo) {
                    lo = min(lo, min(l.a, l.b));
                    hi = max(hi, max(l.a, l.b));
                    fmax = max(fmax, l.factor);
                    count += 1u;
                }
            }
        }
        rf_lo[gi] = float4(lo, fmax);
        rf_hi[gi] = float4(hi, 0.0);
        rf_cs[gi] = uint2(st, count);
        GroupMemoryBarrierWithGroupSync();
        for (uint off = RF_LANES / 2u; off > 0u; off >>= 1u) {
            if (lane < off) {
                float4 l0 = rf_lo[gi], l1 = rf_lo[gi + off];
                float4 h0 = rf_hi[gi], h1 = rf_hi[gi + off];
                uint2  c0 = rf_cs[gi], c1 = rf_cs[gi + off];
                rf_lo[gi] = float4(min(l0.xyz, l1.xyz), max(l0.w, l1.w));
                rf_hi[gi] = float4(max(h0.xyz, h1.xyz), 0.0);
                rf_cs[gi] = uint2(c0.x | c1.x, c0.y + c1.y);
            }
            GroupMemoryBarrierWithGroupSync();
        }
        if (valid && lane == 0u) {
            b.cstate[2u * c] = rf_cs[gi].x;
            if (geo) {
                // The count goes in an integer, not in a float's bits: as a float
                // it is a denormal, which fast math may flush to zero
                b.cstate[2u * c + 1u] = rf_cs[gi].y;
                b.bounds[2u * c + 0u] = rf_lo[gi];
                b.bounds[2u * c + 1u] = rf_hi[gi];
            }
        }
        GroupMemoryBarrierWithGroupSync();                  // the arrays are free for the next eight
    }
}

[shader("compute")]
[numthreads(64, 1, 1)]
void refit(uint3 gid : SV_GroupID, uint gi : SV_GroupIndex) { refit_body<0>(gid, gi); }

[shader("compute")]
[numthreads(64, 1, 1)]
void refit_heavy(uint3 gid : SV_GroupID, uint gi : SV_GroupIndex) { refit_body<1>(gid, gi); }

// The bounding sphere of cluster c of instance inst, in view space, padded by
// the batch's style (scaled with the instance)
struct ClusterSphere {
    float3 c;
    float  r;
    float  pad;                 // of r: the style's (scaled with the instance)
    uint   count;               // visible elements
};

ClusterSphere cluster_sphere(ViewData v, Batch b, uint inst, uint c) {
    ClusterSphere cs;
    float4 lo = b.bounds[2u * c + 0u];
    float4 hi = b.bounds[2u * c + 1u];
    cs.count = b.cstate[2u * c + 1u];
    Mat   T  = load_transform(b.transforms, inst);
    float s  = length(T.c0.xyz);
    float3 ctr = 0.5 * (lo.xyz + hi.xyz);
    float3 ext = 0.5 * (hi.xyz - lo.xyz);
    float pad = b.kind == KIND_SPHERES ? lo.w * b.s0 + b.s1 : b.kind == KIND_SPLINES ? lo.w : lo.w * b.s0;
    cs.c   = mulp(M(v, view), mulp(T, ctr));
    cs.r   = (length(ext) + max(pad, 0.0)) * s;
    cs.pad = max(pad, 0.0) * s;
    return cs;
}

// mode 1: phase 1 over all cluster work. mode 2: phase 2 over the clusters
// phase 1 found occluded. mode 3: the state pass (clusters with a state bit
// in the frustum, occluded or not).
[shader("compute")]
[numthreads(64, 1, 1)]
void cull_clusters(uint3 tid : SV_DispatchThreadID, uint gi : SV_GroupIndex) {
    Args     a = MD_ARGS;
    ViewData v = *a.view;
    Lists    L = *a.lists;
    uint mode = a.mode;
    uint n = mode == 2u ? L.counters[C_COCC] : v.total_cluster_work;
    uint stride = grid_threads(n);
    for (uint i0 = group_start(tid, gi); i0 < n; i0 += stride) {
        uint i = min(i0 + gi, n - 1u);
        bool valid = i0 + gi < n;
        uint slot, inst, c;
        uint2 ent;
        if (mode == 2u) {
            ent  = L.cocc[i];
            slot = (ent.y >> 20) - 1u;
            Batch bb = a.batches[slot];
            uint linear = ent.x;
            inst = linear / max(bb.cluster_count, 1u);
            c    = linear - inst * bb.cluster_count;
        } else {
            slot = find_batch_cluster_work(a.batches, v.batch_count, i);
            Batch bb = a.batches[slot];
            uint local = i - bb.cluster_work_offset;
            inst = local / max(bb.cluster_count, 1u);
            c    = local - inst * bb.cluster_count;
            ent  = uint2(local, (slot + 1u) << 20);
        }
        Batch b = a.batches[slot];
        uint2 rng = cluster_range(b, c);
        ClusterSphere cs = cluster_sphere(v, b, inst, c);
        if (!valid) cs.count = 0u;
        uint hidden = (mode == 1u && valid) ? (rng.y - rng.x) - cs.count : 0u;
        uint out = valid ? OUT_HIDDEN : OUT_NONE;
        // The batch's fade first, over the depths of the cluster's bounding
        // sphere (t grows with depth): gone if its nearest point is past the
        // far band or its farthest short of the near band. Shrinking, its
        // pad shrinks by the largest factor it can have.
        bool faded = false;
        if (cs.count > 0u && mode != 3u && (b.flags & F_FADE) != 0u) {
            float  d  = -cs.c.z;
            float2 tn = fade_ts(b, d - cs.r);
            float2 tf = fade_ts(b, d + cs.r);
            faded = tn.x >= 1.0 || tf.y <= 0.0;
            if ((b.flags & F_FADE_SHRINK) != 0u) cs.r -= cs.pad * (1.0 - (1.0 - tn.x) * tf.y) * 0.999;
        }
        if (cs.count > 0u && faded) {
            out = OUT_LOD;
        } else if (cs.count > 0u) {
            Bounds bd = sphere_bounds(v, cs.c, cs.r);
            if (bd.outside) out = OUT_OUTSIDE;
            else if (mode != 3u && v.occlusion != 0u && hiz_occluded(v, MD_ARGS.hiz, bd, mode == 1u)) out = OUT_OCCLUDED;
            else out = OUT_HW;      // visible: expand
        }
        if (mode == 3u) {
            bool want = valid && cs.count > 0u && b.cstate[2u * c] != 0u;
            uint k = wave_append(L.counters + C_CST, want);
            if (want) L.cst[k] = ent;
            continue;
        }
        uint vbase = mode == 2u ? L.counters[C_CVIS1] : 0u;
        uint iv = wave_append(L.counters + (mode == 2u ? C_CVIS2 : C_CVIS1), out == OUT_HW);
        if (out == OUT_HW) L.cvis[vbase + iv] = ent;
        uint io = wave_append(L.counters + C_COCC, mode == 1u && out == OUT_OCCLUDED);
        if (mode == 1u && out == OUT_OCCLUDED) L.cocc[io] = ent;

        wave_add(L.counters + C_HIDDEN, hidden);
        wave_add(L.counters + C_OUTSIDE, out == OUT_OUTSIDE ? cs.count : 0u);
        wave_add(L.counters + C_LOD, out == OUT_LOD ? cs.count : 0u);
        wave_add(L.counters + C_OCCLUDED, mode == 2u && out == OUT_OCCLUDED ? cs.count : 0u);
        wave_count(L.counters + C_CL_OCC, mode == 2u && out == OUT_OCCLUDED);
        wave_count(L.counters + C_CL_VIS + (mode - 1u), out == OUT_HW);
    }
}

// The fate of one element, in phase 1 or 2 (phase 2: occlusion is final).
// Every lane of a wave calls this together, those without an element with
// valid false: the wave operations inside run converged. (Divergent wave
// operations are legal, and are where drivers have bugs: lavapipe faulted.)
// Texture handles are read from the root where used (MD_ARGS.hiz), never from
// a copy of the argument struct: on Metal, element occlusion through a copied
// struct never fired (the M4 counted no element occluded where lavapipe did).
//
// Spline parts (SPL 1) are not tested for occlusion one by one (their
// clusters are; a part's own test would need a retest list of their own in
// phase 2, for few elements) and have no compute raster.
void element_pass<let SPL : int>(Batch* batches, ViewData v, Lists L, uint phase, bool state_only, bool valid, uint slot, uint inst, uint e) {
    Geo   g;
    g.ok = false; g.faded = false; g.kind = 0u; g.st = 0u; g.a = float3(0.0); g.b = float3(0.0); g.r = 0.0;
    uint ent = 0u;                      // the list entry: global work index
    if (valid) {
        Batch b = batches[slot];
        ent = b.work_offset + inst * b.count + e;
        // the fade: not for the state lists, which do not depend on the
        // camera (the mask's raster applies it)
        if (SPL != 0) g = b.kind == KIND_POLYGONS ? poly_geo(v, b, inst, e, !state_only) : spline_geo(v, b, inst, e, !state_only);
        else          g = element_geo(v, b, inst, e, !state_only);
    }
    uint  out = valid ? OUT_HIDDEN : OUT_NONE;
    if (g.faded) {
        out = OUT_LOD;
    } else if (g.ok) {
        Bounds bd = geo_bounds(v, g);
        float ext = max(bd.pmax.x - bd.pmin.x, bd.pmax.y - bd.pmin.y) - 2.0;     // without the pad
        int4  rp  = reach_px(v, bd);
        if (bd.outside || (!bd.near_clip && ext <= v.small_px && (rp.z < rp.x || rp.w < rp.y))) {
            out = OUT_OUTSIDE;          // off screen, or between pixel centres: no pixel to draw
        } else if (SPL == 0 && !state_only && v.occlusion != 0u && ext > v.occl_min_px && hiz_occluded(v, MD_ARGS.hiz, bd, phase == 1u)) {
            out = OUT_OCCLUDED;
        } else {
            out = (SPL == 0 && v.compute_raster != 0u && !bd.near_clip && ext <= v.small_px) ? OUT_SMALL : OUT_HW;
        }
    }
    if (SPL != 0) {
        bool poly = g.kind == 3u;
        if (state_only) {
            bool want = g.ok && g.st != 0u;
            uint is2 = wave_append(L.counters + C_STS, want && !poly);
            uint is3 = wave_append(L.counters + C_STP, want && poly);
            if (want) { if (poly) L.st3[is3] = ent; else L.st2[is2] = ent; }
            return;
        }
        uint base2 = phase == 2u ? L.counters[C_HWS1] : 0u;
        uint base3 = phase == 2u ? L.counters[C_HWP1] : 0u;
        uint ih2 = wave_append(L.counters + (phase == 2u ? C_HWS2 : C_HWS1), out == OUT_HW && !poly);
        uint ih3 = wave_append(L.counters + (phase == 2u ? C_HWP2 : C_HWP1), out == OUT_HW && poly);
        if (out == OUT_HW) { if (poly) L.hw3[base3 + ih3] = ent; else L.hw2[base2 + ih2] = ent; }
        wave_count(L.counters + C_OUTSIDE, out == OUT_OUTSIDE);
        wave_count(L.counters + C_LOD,     out == OUT_LOD);
        return;
    }
    uint k = g.kind;
    if (state_only) {
        bool want = g.ok && g.st != 0u;
        uint is0 = wave_append(L.counters + C_ST + 0u, want && k == 0u);
        uint is1 = wave_append(L.counters + C_ST + 1u, want && k == 1u);
        if (want) { if (k == 0u) L.st0[is0] = ent; else L.st1[is1] = ent; }
        return;
    }
    uint hw_base0 = phase == 2u ? L.counters[C_HW1 + 0u] : 0u;
    uint hw_base1 = phase == 2u ? L.counters[C_HW1 + 1u] : 0u;
    uint c_hw = phase == 2u ? C_HW2 : C_HW1;
    uint ih0 = wave_append(L.counters + c_hw + 0u, out == OUT_HW && k == 0u);
    uint ih1 = wave_append(L.counters + c_hw + 1u, out == OUT_HW && k == 1u);
    if (out == OUT_HW) { if (k == 0u) L.hw0[hw_base0 + ih0] = ent; else L.hw1[hw_base1 + ih1] = ent; }

    uint sm_base = phase == 2u ? L.counters[C_SMALL1] : 0u;
    uint ism = wave_append(L.counters + (phase == 2u ? C_SMALL2 : C_SMALL1), out == OUT_SMALL);
    if (out == OUT_SMALL) L.small[sm_base + ism] = ent;

    uint ioc = wave_append(L.counters + C_OCC, phase == 1u && out == OUT_OCCLUDED);
    if (phase == 1u && out == OUT_OCCLUDED) L.occ[ioc] = ent;

    // hidden elements were counted with their cluster
    wave_count(L.counters + C_OUTSIDE,  out == OUT_OUTSIDE);
    wave_count(L.counters + C_LOD,      out == OUT_LOD);
    wave_count(L.counters + C_OCCLUDED, phase == 2u && out == OUT_OCCLUDED);
}

// Uniform grid-stride iteration: the group's first item advances for all its
// lanes together, so a lane past the end still takes part (valid false)
uint group_start(uint3 tid, uint gi) { return tid.x - gi; }

// Elements. mode 1: of the clusters visible in phase 1. mode 2: those phase
// 1 found occluded (a thread each). mode 4: of the clusters phase 2 found
// visible. mode 3: the state lists, from the clusters of the state pass.
//
// For the clusters, a group takes 64 of them and its threads walk all their
// elements together (a prefix sum of the sizes, then a binary search per
// element): a cluster of 40 no longer leaves 24 lanes idle, nor one of 70 a
// second round of 58.
groupshared uint4 gc_cl[64];        // slot, instance, first element
groupshared uint  gc_inc[64];       // inclusive prefix of the sizes

void cull_body<let SPL : int>(uint3 tid, uint3 gid, uint gi) {
    Args     a = MD_ARGS;
    ViewData v = *a.view;
    Lists    L = *a.lists;
    uint mode = a.mode;
    if (mode == 2u) {
        if (SPL != 0) return;                               // spline parts and polygons are not in the retest list
        uint n = L.counters[C_OCC];
        uint stride = grid_threads(n);
        for (uint i0 = group_start(tid, gi); i0 < n; i0 += stride) {
            uint i = i0 + gi;
            bool valid = i < n;
            Elem el;
            el.slot = 0u; el.inst = 0u; el.e = 0u;
            if (valid) el = decode_work(a.batches, v.batch_count, L.occ[i]);
            element_pass<0>(a.batches, v, L, 2u, false, valid, el.slot, el.inst, el.e);
        }
        return;
    }
    uint2* src  = mode == 3u ? L.cst : L.cvis;
    uint   base = mode == 4u ? L.counters[C_CVIS1] : 0u;
    uint   n    = mode == 1u ? L.counters[C_CVIS1] : mode == 4u ? L.counters[C_CVIS2] : L.counters[C_CST];
    uint   groups = cluster_groups(n);
    for (uint c0 = gid.x * 64u; c0 < n; c0 += groups * 64u) {
        uint  ci  = c0 + gi;
        uint4 cl  = uint4(0u);
        uint  cnt = 0u;
        if (ci < n) {
            uint2 ent  = src[base + ci];
            uint  slot = (ent.y >> 20) - 1u;
            Batch b    = a.batches[slot];
            uint  inst = ent.x / max(b.cluster_count, 1u);
            uint  c    = ent.x - inst * b.cluster_count;
            uint2 rng  = cluster_range(b, c);
            cl  = uint4(slot, inst, rng.x, 0u);
            cnt = kind_taken<SPL>(b.kind) ? rng.y - rng.x : 0u;     // the other kernel's: none here
        }
        gc_cl[gi]  = cl;
        gc_inc[gi] = cnt;
        GroupMemoryBarrierWithGroupSync();
        for (uint off = 1u; off < 64u; off <<= 1u) {
            uint x = gi >= off ? gc_inc[gi - off] : 0u;
            GroupMemoryBarrierWithGroupSync();
            gc_inc[gi] += x;
            GroupMemoryBarrierWithGroupSync();
        }
        uint total = gc_inc[63];
        for (uint k0 = 0u; k0 < total; k0 += 64u) {         // uniform: element_pass has wave operations
            uint k = k0 + gi;
            bool valid = k < total;
            uint lo = 0u, hi = 63u;                         // the first cluster with gc_inc > k
            while (lo < hi) { uint mid = (lo + hi) >> 1; if (gc_inc[mid] > k) hi = mid; else lo = mid + 1u; }
            uint4 c = gc_cl[lo];
            uint  e = c.z + (k - (lo > 0u ? gc_inc[lo - 1u] : 0u));
            element_pass<SPL>(a.batches, v, L, mode == 4u ? 2u : 1u, mode == 3u, valid, c.x, c.y, e);
        }
        GroupMemoryBarrierWithGroupSync();                  // the arrays are free for the next 64
    }
}

[shader("compute")]
[numthreads(64, 1, 1)]
void cull(uint3 tid : SV_DispatchThreadID, uint3 gid : SV_GroupID, uint gi : SV_GroupIndex) { cull_body<0>(tid, gid, gi); }

[shader("compute")]
[numthreads(64, 1, 1)]
void cull_heavy(uint3 tid : SV_DispatchThreadID, uint3 gid : SV_GroupID, uint gi : SV_GroupIndex) { cull_body<1>(tid, gid, gi); }

// The indirect commands and grids made of the counts. One thread.
//   11: after the phase 1 cluster cull: element grid, phase 2 cluster grid
//   12: after the phase 1 element cull: draws, small grid, element retest grid
//   21: after the phase 2 cluster cull: element grid
//   22: after the phase 2 element culls: draws, small grid
//   31: after the state cluster pass: its element grid
//   32: after the state element pass: the state draws
[shader("compute")]
[numthreads(1, 1, 1)]
void finalize(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    uint* c = a.lists->counters;
    uint mode = a.mode;
    if (mode == 11u || mode == 21u || mode == 31u) {
        uint n = mode == 11u ? c[C_CVIS1] : mode == 21u ? c[C_CVIS2] : c[C_CST];
        uint* g = c + (mode == 31u ? C_GRID_CST : C_GRID_CVIS + (mode == 21u ? 4u : 0u));
        g[0] = cluster_groups(n); g[1] = 1u; g[2] = 1u;
        if (mode == 11u) {
            uint no = c[C_COCC];
            uint* go = c + C_GRID_COCC;
            go[0] = min((no + 63u) / 64u, 65535u); go[1] = 1u; go[2] = 1u;
        }
        return;
    }
    if (mode == 32u) {
        uint* sc = c + C_STCMD;
        sc[0] = 6u; sc[1] = c[C_ST]; sc[2] = 0u; sc[3] = 0u;
        uint* si = c + C_STICMD;                // capsules: indexed, 36 indices over 8 corners
        si[0] = 36u; si[1] = c[C_ST + 1u]; si[2] = 0u; si[3] = 0u; si[4] = 0u;
        uint* ss = c + C_STICMDS;               // spline parts: the same boxes
        ss[0] = 36u; ss[1] = c[C_STS]; ss[2] = 0u; ss[3] = 0u; ss[4] = 0u;
        uint* sp = c + C_STICMDP;               // polygons: boxes too
        sp[0] = 36u; sp[1] = c[C_STP]; sp[2] = 0u; sp[3] = 0u; sp[4] = 0u;
        return;
    }
    uint phase = mode == 12u ? 1u : 2u;
    uint* cmd = c + C_CMD + (phase - 1u) * 2u * 4u;
    cmd[0] = 6u;
    cmd[1] = phase == 1u ? c[C_HW1] : c[C_HW2];
    cmd[2] = 0u;
    cmd[3] = phase == 1u ? 0u : c[C_HW1];
    uint* icmd = c + C_ICMD + (phase - 1u) * 5u;    // capsules: indexed, 36 indices over 8 corners
    icmd[0] = 36u;
    icmd[1] = phase == 1u ? c[C_HW1 + 1u] : c[C_HW2 + 1u];
    icmd[2] = 0u;
    icmd[3] = 0u;
    icmd[4] = phase == 1u ? 0u : c[C_HW1 + 1u];
    uint* scmd = c + C_ICMDS + (phase - 1u) * 5u;   // spline parts: the same boxes
    scmd[0] = 36u;
    scmd[1] = phase == 1u ? c[C_HWS1] : c[C_HWS2];
    scmd[2] = 0u;
    scmd[3] = 0u;
    scmd[4] = phase == 1u ? 0u : c[C_HWS1];
    uint* pcmd = c + C_ICMDP + (phase - 1u) * 5u;   // polygons: boxes too
    pcmd[0] = 36u;
    pcmd[1] = phase == 1u ? c[C_HWP1] : c[C_HWP2];
    pcmd[2] = 0u;
    pcmd[3] = 0u;
    pcmd[4] = phase == 1u ? 0u : c[C_HWP1];
    uint ns = phase == 1u ? c[C_SMALL1] : c[C_SMALL2];
    uint* gs = c + C_GRID_SMALL + (phase - 1u) * 4u;
    gs[0] = min((ns + 63u) / 64u, 65535u); gs[1] = 1u; gs[2] = 1u;
    if (phase == 1u) {
        uint no = c[C_OCC];
        uint* go = c + C_GRID_OCC;
        go[0] = min((no + 63u) / 64u, 65535u); go[1] = 1u; go[2] = 1u;
    }
}

// ---------------------------------------------------------------------------
// Compute raster of small primitives
// ---------------------------------------------------------------------------

// 64-bit atomic max without a return value: Metal has only that form for
// 64-bit integers (atomic_max_explicit on atomic_ulong, Apple8+), and Slang's
// Atomic<uint64_t>.max() lowers to the fetching one there.
void atomic_max_u64(uint64_t* p, uint64_t v) {
    __target_switch {
    case metal:
        __intrinsic_asm "atomic_max_explicit((device atomic_ulong*)($0), $1, memory_order_relaxed)";
    default:
        ((Atomic<uint64_t>*)p)->max(v);
    }
}

// Every pixel a small primitive can cover is tested against the exact shape,
// with the same ray and intersection as the fragment stages, and the nearest
// wins by a 64-bit max of (depth, index into the small list). Depth bits
// compare as unsigned because depths are non-negative.
//
// A group takes 64 primitives and shares out all their pixels among its 64
// threads (a prefix sum of the pixel counts, then a binary search per pixel),
// so a lane with a large primitive does not hold up 63 with small ones. The
// pixels are those whose centre the jittered shape can reach (reach_px: the
// exact rectangle widened by this raster's jitter). A pixel already holding a
// nearer depth skips the atomic (the read is a plain one; a stale value only
// costs an atomic that changes nothing).
groupshared float4 gs_a[64];        // a, radius
groupshared float4 gs_b[64];        // b, kind
groupshared uint4  gs_box[64];      // x0, y0, width, index into the small list
groupshared uint   gs_inc[64];      // inclusive prefix of the pixel counts

[shader("compute")]
[numthreads(64, 1, 1)]
void raster_small(uint3 tid : SV_DispatchThreadID, uint gi : SV_GroupIndex) {
    Args     a = MD_ARGS;
    ViewData v = *a.view;
    Lists    L = *a.lists;
    uint phase = a.mode;
    uint n     = phase == 1u ? L.counters[C_SMALL1] : L.counters[C_SMALL2];
    uint base  = phase == 1u ? 0u : L.counters[C_SMALL1];
    uint stride = grid_threads(n);
    int  W = int(v.viewport.x), H = int(v.viewport.y);
    uint* depth_words = (uint*)L.vis64;
    for (uint i0 = group_start(tid, gi); i0 < n; i0 += stride) {
        uint i   = i0 + gi;
        uint cnt = 0u;
        uint idx = base + i;
        float4 ga = float4(0.0), gb = float4(0.0);
        uint4  box = uint4(0u, 0u, 1u, idx);
        if (i < n) {
            Elem  el = decode_work(a.batches, v.batch_count, L.small[idx]);
            Batch b  = a.batches[el.slot];
            Geo   g  = element_geo(v, b, el.inst, el.e, true);
            if (g.ok) {
                Bounds bd = geo_bounds(v, g);
                int4 rp = reach_px(v, bd);
                int2 p0 = max(rp.xy, int2(0));
                int2 p1 = min(rp.zw, int2(W - 1, H - 1));
                if (p1.x >= p0.x && p1.y >= p0.y) {
                    uint w = uint(p1.x - p0.x + 1), h = uint(p1.y - p0.y + 1);
                    cnt = w * h;
                    box = uint4(uint(p0.x), uint(p0.y), w, idx);
                    ga  = float4(g.a, g.r);
                    gb  = float4(g.b, float(g.kind));
                }
            }
        }
        gs_a[gi]   = ga;
        gs_b[gi]   = gb;
        gs_box[gi] = box;
        gs_inc[gi] = cnt;
        GroupMemoryBarrierWithGroupSync();
        for (uint off = 1u; off < 64u; off <<= 1u) {
            uint x = gi >= off ? gs_inc[gi - off] : 0u;
            GroupMemoryBarrierWithGroupSync();
            gs_inc[gi] += x;
            GroupMemoryBarrierWithGroupSync();
        }
        uint total = gs_inc[63];
        for (uint k = gi; k < total; k += 64u) {
            uint lo = 0u, hi = 63u;                         // the first e with gs_inc[e] > k
            while (lo < hi) { uint mid = (lo + hi) >> 1; if (gs_inc[mid] > k) hi = mid; else lo = mid + 1u; }
            uint  e     = lo;
            uint  local = k - (e > 0u ? gs_inc[e - 1u] : 0u);
            uint4 bx    = gs_box[e];
            uint  x     = bx.x + local % bx.z;
            uint  y     = bx.y + local / bx.z;
            float4 ea   = gs_a[e];
            float4 eb   = gs_b[e];
            Ray   r = pixel_ray(v, float2(float(x), float(y)) + 0.5);
            float t;
            if (eb.w < 0.5) t = sphere_hit(r, ea.xyz, ea.w, false);
            else            t = capsule_hit(r, ea.xyz, eb.xyz, ea.w);
            if (t < 0.0) continue;
            float d = clip_depth(v, r.o + r.d * t);
            if (eb.w >= 0.5) d = capsule_depth(d);
            uint pix = y * uint(W) + x;
            if (depth_words[2u * pix + 1u] > asuint(d)) continue;      // nearer already there
            atomic_max_u64(L.vis64 + pix, (uint64_t(asuint(d)) << 32) | uint64_t(bx.w));
        }
        GroupMemoryBarrierWithGroupSync();                  // the arrays are free for the next 64
    }
}

// The compute raster's result into depth and ids, so that it occludes and is
// occluded like any other raster. A fullscreen triangle.
struct FullV {
    float4 pos : SV_Position;
};

[shader("vertex")]
FullV merge_vs(uint vid : SV_VertexID) {
    FullV o;
    float2 p = float2(float((vid << 1) & 2u), float(vid & 2u));
    o.pos = float4(p * 2.0 - 1.0, 0.0, 1.0);
    if (MD_ARGS.width == 0u) o.pos = float4(0.0, 0.0, -2.0, 1.0);  // reads the arguments: both stages must declare the same struct
    return o;
}

struct MergeOut {
    uint2 id    : SV_Target0;
    float depth : SV_Depth;
};

[shader("fragment")]
MergeOut merge_fs(FullV i) {
    Args a = MD_ARGS;
    uint2 px = uint2(i.pos.xy);
    uint64_t val = a.lists->vis64[px.y * a.width + px.x];
    if (val == 0ull) discard;
    MergeOut o;
    o.depth = asfloat(uint(val >> 32));
    o.id    = elem_id(a.batches, decode_work(a.batches, a.view->batch_count, a.lists->small[uint(val & 0xFFFFFFFFull)]));
    return o;
}

// ---------------------------------------------------------------------------
// Depth pyramid
// ---------------------------------------------------------------------------

// One level from the one before (or from the depth target): the farthest
// (min) depth of each 2 x 2, edges clamped, so every texel covers exactly its
// pixels.
[shader("compute")]
[numthreads(8, 8, 1)]
void hiz_reduce(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    if (tid.x >= a.dst_w || tid.y >= a.dst_h) return;
    int2 s  = int2(tid.xy) * 2;
    int2 mx = int2(int(a.src_w) - 1, int(a.src_h) - 1);
    float m = a.src.Load(int3(min(s, mx), int(a.src_mip)));
    m = min(m, a.src.Load(int3(min(s + int2(1, 0), mx), int(a.src_mip))));
    m = min(m, a.src.Load(int3(min(s + int2(0, 1), mx), int(a.src_mip))));
    m = min(m, a.src.Load(int3(min(s + int2(1, 1), mx), int(a.src_mip))));
    a.dst[tid.xy] = m;
}

// ---------------------------------------------------------------------------
// Hardware raster, pulling from a list
// ---------------------------------------------------------------------------

static const float4 CULLED = float4(0.0, 0.0, -2.0, 1.0);   // outside the depth range: clipped

uint list_entry(Args a, uint i) {
    Lists L = *a.lists;
    switch (a.list) {
    case 0u: return L.hw0[i];
    case 1u: return L.hw1[i];
    case 2u: return L.st0[i];
    case 3u: return L.st1[i];
    case 4u: return L.hw2[i];
    case 5u: return L.st2[i];
    case 6u: return L.hw3[i];
    default: return L.st3[i];
    }
}

struct SphereV {
    float4 pos : SV_Position;
    nointerpolation float4 sph   : SPH;     // view-space centre, radius
    nointerpolation uint2  idv   : IDV;
    nointerpolation uint   state : STATE;
};

[shader("vertex")]
SphereV sphere_vs(uint vid : SV_VertexID, uint iid : SV_InstanceID,
                  uint base_vertex : SV_StartVertexLocation, uint base_instance : SV_StartInstanceLocation) {
    Args     a = MD_ARGS;
    ViewData v = *a.view;
    SphereV  o;
    o.pos   = CULLED;
    o.sph   = float4(0.0);
    o.idv   = uint2(0u);
    o.state = 0u;

    uint  corner = md_draw_vertex(vid, base_vertex);
    Elem  el     = decode_work(a.batches, v.batch_count, list_entry(a, base_instance + md_draw_instance(iid, base_instance)));
    uint2 ent    = elem_id(a.batches, el);
    Geo   g      = element_geo(v, a.batches[el.slot], el.inst, el.e, true);   // the state lists' fade here
    if (!g.ok) return o;
    float3 c   = g.a;
    float  rad = g.r;

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
    // Proxies that only touch the silhouette lose pixels whose centres lie on
    // it to rounding (the compute raster, which tests every pixel of a padded
    // rectangle, showed it): a little larger, still in front of the surface.
    qr *= 1.002;
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
    o.idv   = ent;
    o.state = g.st;
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

// The box around a capsule: drawn indexed (md_render.c holds the 36 indices,
// 12 triangles counter-clockwise seen from outside), so the vertex id is the
// corner c = (c & 1, c >> 1 & 1, c >> 2 & 1) in the capsule's frame
// (u, v, axis), and each corner is shaded once.

[shader("vertex")]
CapsuleV capsule_vs(uint vid : SV_VertexID, uint iid : SV_InstanceID,
                    uint base_vertex : SV_StartVertexLocation, uint base_instance : SV_StartInstanceLocation) {
    Args     a = MD_ARGS;
    ViewData v = *a.view;
    CapsuleV o;
    o.pos   = CULLED;
    o.pa    = float4(0.0);
    o.pb    = float4(0.0);
    o.idv   = uint2(0u);
    o.state = 0u;

    uint  vk  = md_draw_vertex(vid, base_vertex);
    Elem  el  = decode_work(a.batches, v.batch_count, list_entry(a, base_instance + md_draw_instance(iid, base_instance)));
    uint2 ent = elem_id(a.batches, el);
    Geo   g   = element_geo(v, a.batches[el.slot], el.inst, el.e, true);   // the state lists' fade here
    if (!g.ok) return o;
    float4 corner = capsule_box_corner(v, g, vk);
    o.pos   = corner;
    o.pa    = float4(g.a, g.r);
    o.pb    = float4(g.b, 0.0);
    o.idv   = ent;
    o.state = g.st;
    return o;
}

// Corner c of the box around the capsule g, in clip space
float4 capsule_box_corner(ViewData v, Geo g, uint c) {
    float3 pa  = g.a;
    float3 pb  = g.b;
    float  rad = g.r;
    float3 ax  = pb - pa;
    float  len = length(ax);
    float3 z   = len > 1e-20 ? ax / len : float3(0.0, 0.0, 1.0);
    float3 hlp = abs(z.x) < 0.9 ? float3(1.0, 0.0, 0.0) : float3(0.0, 1.0, 0.0);
    float3 x   = normalize(cross(z, hlp));
    float3 y   = cross(z, x);               // (x, y, z) right-handed

    // Built in view space from the transformed ends, in a right-handed frame:
    // outward faces stay counter-clockwise whatever the transforms mirror.
    float  br = rad * 1.002;                // see sphere_vs: not just touching
    float3 p = pa
             + x * (((c & 1u) != 0u) ? br : -br)
             + y * (((c & 2u) != 0u) ? br : -br)
             + z * (((c & 4u) != 0u) ? len + br : -br);
    return mul4(M(v, proj), float4(p, 1.0));
}

[shader("fragment")]
VisOut capsule_fs(CapsuleV i) {
    ViewData v = *MD_ARGS.view;
    Ray   r = pixel_ray(v, i.pos.xy);
    float t = capsule_hit(r, i.pa.xyz, i.pb.xyz, i.pa.w);
    if (t < 0.0) discard;
    VisOut o;
    o.id    = i.idv;
    o.depth = capsule_depth(clip_depth(v, r.o + r.d * t));
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
    o.depth = capsule_depth(clip_depth(v, r.o + r.d * t));
    return o;
}

// Spline parts: the box around the part's bounding capsule (element_geo's
// for capsules, spline_geo's here), a fragment casting its ray against the
// part rebuilt from its control points, cut into the pieces the vertex
// stage counted (as eval_hit counts them)
struct SplineV {
    float4 pos : SV_Position;
    nointerpolation uint4  pe    : PE;      // pieces, instance, element
    nointerpolation uint2  idv   : IDV;
    nointerpolation uint   state : STATE;
};

[shader("vertex")]
SplineV spline_vs(uint vid : SV_VertexID, uint iid : SV_InstanceID,
                  uint base_vertex : SV_StartVertexLocation, uint base_instance : SV_StartInstanceLocation) {
    Args     a = MD_ARGS;
    ViewData v = *a.view;
    SplineV  o;
    o.pos   = CULLED;
    o.pe    = uint4(0u);
    o.idv   = uint2(0u);
    o.state = 0u;
    uint  vk  = md_draw_vertex(vid, base_vertex);
    Elem  el  = decode_work(a.batches, v.batch_count, list_entry(a, base_instance + md_draw_instance(iid, base_instance)));
    Batch b   = a.batches[el.slot];
    Geo   g   = spline_geo(v, b, el.inst, el.e, true);     // the state lists' fade here
    if (!g.ok) return o;
    o.pos   = capsule_box_corner(v, g, vk);
    o.pe    = uint4(spline_element_pieces(v, b, el.inst, el.e, g), el.inst, el.e, el.slot);
    o.idv   = elem_id(a.batches, el);
    o.state = g.st;
    return o;
}

float spline_fragment_hit(ViewData v, Ray r, SplineV i) {
    return spline_hit(r, spline_part(v, MD_ARGS.batches[i.pe.w], i.pe.y, i.pe.z, true), i.pe.x).t;
}

[shader("fragment")]
VisOut spline_fs(SplineV i) {
    ViewData v = *MD_ARGS.view;
    Ray   r = pixel_ray(v, i.pos.xy);
    float t = spline_fragment_hit(v, r, i);
    if (t < 0.0) discard;
    VisOut o;
    o.id    = i.idv;
    o.depth = capsule_depth(clip_depth(v, r.o + r.d * t));
    return o;
}

[shader("fragment")]
MaskOut spline_mask_fs(SplineV i) {
    ViewData v = *MD_ARGS.view;
    Ray   r = pixel_ray(v, i.pos.xy);
    float t = spline_fragment_hit(v, r, i);
    if (t < 0.0) discard;
    MaskOut o;
    o.m     = float4(0.0, (i.state & 1u) != 0u ? 1.0 : 0.0, 0.0, (i.state & 2u) != 0u ? 1.0 : 0.0);
    o.depth = capsule_depth(clip_depth(v, r.o + r.d * t));
    return o;
}

// Polygons: the box around the polygon (poly_box_corner), a fragment casting
// its ray against the polygon rebuilt from its corners. SplineV's layout,
// pe.x unused.
[shader("vertex")]
SplineV polygon_vs(uint vid : SV_VertexID, uint iid : SV_InstanceID,
                   uint base_vertex : SV_StartVertexLocation, uint base_instance : SV_StartInstanceLocation) {
    Args     a = MD_ARGS;
    ViewData v = *a.view;
    SplineV  o;
    o.pos   = CULLED;
    o.pe    = uint4(0u);
    o.idv   = uint2(0u);
    o.state = 0u;
    uint  vk  = md_draw_vertex(vid, base_vertex);
    Elem  el  = decode_work(a.batches, v.batch_count, list_entry(a, base_instance + md_draw_instance(iid, base_instance)));
    Batch b   = a.batches[el.slot];
    PolyPart pp;
    Geo   g   = poly_geo_part(v, b, el.inst, el.e, true, pp);     // the state lists' fade here
    if (!g.ok) return o;
    o.pos   = poly_box_corner(v, pp, vk);
    o.pe    = uint4(0u, el.inst, el.e, el.slot);
    o.idv   = elem_id(a.batches, el);
    o.state = g.st;
    return o;
}

float poly_fragment_hit(ViewData v, Ray r, SplineV i) {
    return poly_hit(r, poly_part(v, MD_ARGS.batches[i.pe.w], i.pe.y, i.pe.z, true)).t;
}

[shader("fragment")]
VisOut polygon_fs(SplineV i) {
    ViewData v = *MD_ARGS.view;
    Ray   r = pixel_ray(v, i.pos.xy);
    float t = poly_fragment_hit(v, r, i);
    if (t < 0.0) discard;
    VisOut o;
    o.id    = i.idv;
    o.depth = capsule_depth(clip_depth(v, r.o + r.d * t));
    return o;
}

[shader("fragment")]
MaskOut polygon_mask_fs(SplineV i) {
    ViewData v = *MD_ARGS.view;
    Ray   r = pixel_ray(v, i.pos.xy);
    float t = poly_fragment_hit(v, r, i);
    if (t < 0.0) discard;
    MaskOut o;
    o.m     = float4(0.0, (i.state & 1u) != 0u ? 1.0 : 0.0, 0.0, (i.state & 2u) != 0u ? 1.0 : 0.0);
    o.depth = capsule_depth(clip_depth(v, r.o + r.d * t));
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
            uint w = b.palette_index[gated(b, F_PALETTE_IDX, ai) >> 1];
            idx = ((ai & 1u) != 0u) ? (w >> 16) : (w & 0xFFFFu);
        }
        if (idx >= b.palette_count) idx = 0u;
        c = b.palette_count > 0u ? b.palette[idx] : 0xFFFFFFFFu;     // idx < palette_count, else 0
    } else {
        c = ((b.flags & F_COLOR) != 0u) ? b.color[gated(b, F_COLOR, ai)] : 0xFFFFFFFFu;
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

// SPL 0: spheres and capsules; 1: splines and polygons (the other's pixels
// come back invalid)
Hit eval_hit<let SPL : int>(ViewData v, Batch* batches, uint2 idv, float2 pix) {
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
    h.id       = ((b.flags & F_ID) != 0u) ? b.id[gated(b, F_ID, h.element)] : h.element;
    h.attr     = 0xFFFFFFFFu;
    if (b.kind == KIND_EXTERNAL) return h;

    Mat   T  = load_transform(b.transforms, h.instance);
    Mat   Tp = load_transform(b.transforms_prev, h.instance);
    float s  = length(T.c0.xyz);
    Ray   r  = pixel_ray(v, pix);

    if (SPL == 0 && b.kind == KIND_SPHERES) {
        uint   pi  = point_of(b, h.element, 0u, 1u);
        float3 lc  = load_pos(b, pi);
        float3 cw  = mulp(T, lc);
        float3 cwp = mulp(Tp, lc - load_vel(b, pi) * v.dt);
        float3 c   = mulp(M(v, view), cw);
        float  rad = sphere_radius(b, pi) * s;
        if ((b.flags & F_FADE) != 0u && (b.flags & F_FADE_SHRINK) != 0u) rad *= fade_size(fade_ts(b, -c.z));   // as element_geo
        float  t   = sphere_hit(r, c, rad, true);
        h.p_view   = r.o + r.d * t;
        h.n_view   = normalize(h.p_view - c);
        h.p_world  = mulp(M(v, inv_view), h.p_view);
        h.p_world_prev = h.p_world - (cw - cwp);
        h.attr     = attr_of(b, pi);
        h.color    = attr_color(b, h.attr);
        h.valid    = true;
    } else if (SPL == 0 && b.kind == KIND_CAPSULES) {
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
        if ((b.flags & F_FADE) != 0u && (b.flags & F_FADE_SHRINK) != 0u) rad *= fade_size(fade_ts(b, -0.5 * (pa.z + pb.z)));
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
    } else if (SPL != 0 && b.kind == KIND_SPLINES) {
        // the part again, cut as the raster cut it (its bounds give the pieces;
        // not through spline_geo, which crashed lavapipe's pick kernel)
        SplinePart sp = spline_part(v, b, h.instance, h.element, true);
        Geo        g;
        g.ok = true; g.faded = false; g.kind = 1u; g.st = 0u;
        spline_bounds(sp, g.a, g.b, g.r);                // as element_geo: the same pieces
        SplineHit  sh = spline_hit(r, sp, spline_pieces(v, sp, g));
        if (sh.t < 0.0) {                               // a rounding miss at the silhouette: the nearest point to its control point
            sh.t = max(dot(sp.p2 - r.o, r.d), 0.0);
            sh.n = -r.d;
        }
        uint   pi  = b.index[h.element] & SPLINE_POINT;
        float3 lc  = load_pos(b, pi);
        float3 cw  = mulp(T, lc);
        float3 cwp = mulp(Tp, lc - load_vel(b, pi) * v.dt);
        h.p_view   = r.o + r.d * sh.t;
        h.n_view   = sh.n;
        h.p_world  = mulp(M(v, inv_view), h.p_view);
        h.p_world_prev = h.p_world - (cw - cwp);
        h.attr     = attr_of(b, pi);
        h.color    = attr_color(b, h.attr);
        h.valid    = true;
    } else if (SPL != 0 && b.kind == KIND_POLYGONS) {
        PolyPart  pp = poly_part(v, b, h.instance, h.element, true);
        SplineHit ph = poly_hit(r, pp);
        if (ph.t < 0.0) {                               // a rounding miss at the silhouette: the nearest point to its centroid
            ph.t = max(dot(pp.c - r.o, r.d), 0.0);
            ph.n = -r.d;
        }
        float3 dm = float3(0.0);                        // its motion: its corners' mean
        for (uint k = 0u; k < pp.n; ++k) {
            uint   pk = point_of(b, h.element, k, pp.n);
            float3 lk = load_pos(b, pk);
            dm += mulp(T, lk) - mulp(Tp, lk - load_vel(b, pk) * v.dt);
        }
        h.p_view   = r.o + r.d * ph.t;
        h.n_view   = ph.n;
        h.p_world  = mulp(M(v, inv_view), h.p_view);
        h.p_world_prev = h.p_world - dm / float(pp.n);
        uint  kn = 0u;                                  // its attribute: the corner nearest the hit
        float dn = 3.0e38;
        [ForceUnroll] for (uint k = 0u; k < POLY_MAX; ++k) {
            if (k < pp.n) { float d = length(pp.q[k] - h.p_view); if (d < dn) { dn = d; kn = k; } }
        }
        h.attr     = attr_of(b, point_of(b, h.element, kn, pp.n));
        h.color    = attr_color(b, attr_of(b, point_of(b, h.element, 0u, pp.n)));   // its colour: the first's
        h.valid    = true;
    }
    return h;
}

float2 to_pixels(ViewData v, float4 clip) {
    float2 ndc = clip.xy / clip.w;
    return float2((ndc.x * 0.5 + 0.5) * v.viewport.x, (0.5 - ndc.y * 0.5) * v.viewport.y);
}

// Each pixel by the twin that takes its kind: the main one clears the empty
// ones, the heavy twin (run after it, when the frame has splines or
// polygons) passes over all but theirs
void resolve_body<let SPL : int>(uint3 tid) {
    Args a = MD_ARGS;
    if (tid.x >= a.width || tid.y >= a.height) return;
    int2  px  = int2(tid.xy);
    uint2 idv = a.id.Load(int3(px, 0));
    if ((idv.y >> 20) == 0u) {
        if (SPL != 0) return;
        a.color[px]  = float4(0.0);
        a.normal[px] = float2(0.0);
        a.motion[px] = float2(0.0);
        return;
    }
    if (!kind_taken<SPL>(a.batches[(idv.y >> 20) - 1u].kind)) return;
    ViewData v = *a.view;
    Hit h = eval_hit<SPL>(v, a.batches, idv, float2(tid.xy) + 0.5);
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

[shader("compute")]
[numthreads(8, 8, 1)]
void resolve(uint3 tid : SV_DispatchThreadID) { resolve_body<0>(tid); }

[shader("compute")]
[numthreads(8, 8, 1)]
void resolve_heavy(uint3 tid : SV_DispatchThreadID) { resolve_body<1>(tid); }

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
        if (b.kind == KIND_SPHERES || b.kind == KIND_CAPSULES || b.kind == KIND_POLYGONS || b.kind == KIND_SPLINES) {
            uint64_t linear = (uint64_t(idv.y & 0x3FFFFu) << 32) | uint64_t(idv.x);
            uint cnt = max(b.count, 1u);
            uint e   = uint(linear % uint64_t(cnt));
            uint st;
            if (b.kind == KIND_SPHERES) {
                st = attr_state(b, attr_of(b, point_of(b, e, 0u, 1u)));
            } else if (b.kind == KIND_SPLINES) {
                st = attr_state(b, attr_of(b, b.index[e] & SPLINE_POINT));
            } else if (b.kind == KIND_POLYGONS) {
                poly_visible(b, e, st);
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
    Args  a = MD_ARGS;
    uint* o = a.pick_out;
    for (uint i = 0u; i < 12u; ++i) o[i] = 0u;
    uint px = a.px, py = a.py;
    if (px >= a.width || py >= a.height) return;
    uint2 idv = a.id.Load(int3(int(px), int(py), 0));
    o[0] = idv.x;
    o[1] = idv.y;
    if ((idv.y >> 20) == 0u) return;
    ViewData v   = *a.view;
    float2   pix = float2(float(px), float(py)) + 0.5;
    Hit h;
    if (heavy_kind(a.batches[(idv.y >> 20) - 1u].kind)) h = eval_hit<1>(v, a.batches, idv, pix);
    else                                                h = eval_hit<0>(v, a.batches, idv, pix);
    o[2] = h.id;
    o[3] = h.attr;
    if (h.valid) {
        o[4] = asuint(h.p_world.x);
        o[5] = asuint(h.p_world.y);
        o[6] = asuint(h.p_world.z);
        o[7] = asuint(clip_depth(v, h.p_view));
    }
}
