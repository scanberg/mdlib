// md_compose.slang -- the kernels of md_compose.c: the ambient occlusion of
// the opaque surface (viamd's ssao/ssao.frag, blur.frag, upsample.frag,
// depth_downsample.frag), and the compose: the opaque surface shaded under
// it (viamd's compose_deferred.frag, both models), the transparent layers
// shaded and composited over it front to back.
//
// Pixels are x right, y down, their centres at + 0.5 (md_render's). View
// space looks down -z; a linear depth is -z. Layout rules of md_gpu.slang: no
// 3-vectors, 4-vectors at 16-byte offsets, pointers and handles at 8.
// Mirrored in md_compose.c.

#include "../md_gpu.slang"

struct Args {
    float4 ip0, ip1, ip2, ip3;  // the inverse projection's columns (the layout check models no arrays)
    float4 light;               // xyz: toward the key light (unit), w: its strength
    float4 dir_radiance;        // xyz, w: F0
    float4 sky;                 // xyz, w: roughness
    float4 ground;              // xyz, w: basic ambient
    float4 misc;                // x: basic specular, y: 1 basic, z: 1 linear colours, w: 1 orthographic
    float4 pa;                  // the projection: P00, P11, P20, P21 (Pcr: column c, row r)
    float4 pb;                  // P30, P31, P22, P32
    float4 pc;                  // P23, P33, a pixel's size in view space at linear depth 1 (or at any:
                                // orthographic), the occlusion's intensity
    uint*   layers;             // md_render_targets_t.layers: colour, normal, depth planes
    float*  lin;                // the linear depth's levels 0..AO_LEVELS - 1, one after another
    float2* ao;                 // two planes of the half resolution (level 1): (visibility, linear depth)
    uint   count;               // layers a pixel (K)
    uint   stride;              // words of a layer (pixels)
    uint   width;
    uint   height;
    uint   level;               // ao_linear: the level made
    uint   dir;                 // ao_blur: 0 along x, 1 along y
    uint   src;                 // ao_blur: the plane read (the other is written)
    uint   use_ao;              // compose: plane 0 holds the occlusion
    DescriptorHandle<RWTexture2D<float4>> target;
    DescriptorHandle<Texture2D<float>>    depth;    // md_render's: reverse Z, 0 where nothing
    DescriptorHandle<Texture2D<float4>>   color;    // sRGB-encoded base colour, a: covered
    DescriptorHandle<Texture2D<float2>>   normal;   // view space, octahedral
    float4* vols;               // the volumes' records, VOL_STRIDE float4 each (md_compose.c's vol_record)
    uint   vol_count;
    uint   frame;               // seeds the emission's jitter
    DescriptorHandle<Texture3D<float>>    vf0, vf1, vf2, vf3;   // the volumes' fields
    DescriptorHandle<Texture2D<float4>>   vt0, vt1, vt2, vt3;   // and transfer functions
    DescriptorHandle<SamplerState>        trilinear;            // clamp to edge
    DescriptorHandle<Texture3D<float4>>   vc0, vc1, vc2, vc3;   // colour volumes
    DescriptorHandle<Texture3D<float>>    vg0, vg1, vg2, vg3;   // colour fields
    DescriptorHandle<Texture2D<float4>>   vm0, vm1, vm2, vm3;   // their colour maps
    DescriptorHandle<Texture3D<float4>>   vn0, vn1, vn2, vn3;   // normals
    float4* slices;             // the slices' records, SLICE_STRIDE float4 each (md_compose.c's)
    uint   slice_count;
    uint   _sp;
    DescriptorHandle<Texture3D<float>>    sf0, sf1, sf2, sf3, sf4, sf5, sf6, sf7;   // the slices' fields
    DescriptorHandle<Texture2D<float4>>   st0, st1, st2, st3, st4, st5, st6, st7;   // and transfer functions
};
MD_SHADER_ARGS(Args);

static const float PI          = 3.1415926535;
static const float ONE_OVER_PI = 1.0 / 3.1415926535;

// ---- Views ------------------------------------------------------------------

static const uint  AO_LEVELS = 6u;      // 0 .. AO_MAX_MIP
static const float Z_BG      = 1.0e9;   // the linear depth where nothing was drawn
static const float Z_MAX     = 0.5e9;   // at or beyond: nothing

uint2 level_size(Args a, uint l) { return uint2(max(a.width >> l, 1u), max(a.height >> l, 1u)); }

uint level_offset(Args a, uint l) {
    uint o = 0u;
    for (uint k = 0u; k < l; ++k) { uint2 s = level_size(a, k); o += s.x * s.y; }
    return o;
}

float lin_at(Args a, uint l, int2 q) {
    uint2 s = level_size(a, l);
    uint2 c = uint2(clamp(q, int2(0), int2(s) - 1));
    return a.lin[level_offset(a, l) + c.y * s.x + c.x];
}

float2 ndc_of(Args a, float2 px) {
    return float2(px.x / float(a.width) * 2.0 - 1.0, 1.0 - px.y / float(a.height) * 2.0);
}

// The view-space point at a pixel (its coordinates, not its centre's index)
// and linear depth z
float3 view_at(Args a, float2 px, float z) {
    float2 n  = ndc_of(a, px);
    float  zv = -z;
    float  w  = a.pc.x * zv + a.pc.y;
    return float3((n.x * w - a.pa.z * zv - a.pb.x) / a.pa.x, (n.y * w - a.pa.w * zv - a.pb.y) / a.pa.y, zv);
}

float linear_depth(Args a, float d) {
    if (d == 0.0) return Z_BG;
    float zv = (a.pb.w - d * a.pc.y) / (d * a.pc.x - a.pb.z);
    return -zv;
}

// A pixel's size in view space at linear depth z
float px_world(Args a, float z) { return a.misc.w > 0.5 ? a.pc.z : a.pc.z * z; }

float3 oct_decode(float2 e) {
    float3 n = float3(e, 1.0 - abs(e.x) - abs(e.y));
    if (n.z < 0.0) {
        float2 s = select(n.xy >= 0.0, float2(1.0), float2(-1.0));
        n.xy = (1.0 - abs(n.yx)) * s;
    }
    return normalize(n);
}

float3 normal_at(Args a, int2 q) { return oct_decode(a.normal.Load(int3(q, 0))); }

// A half float's bits as a float (finite values; f16tof32 declares Int16,
// which not every device has): its bits moved into a float's, rebiased
float half_to_float(uint h) {
    float f = asfloat((h & 0x7FFFu) << 13) * 5.192296858534828e33;   // 2^112
    return (h & 0x8000u) != 0u ? -f : f;
}

// The full resolution pixel a texel of level 1 holds (the rotated grid's pick)
int2 pick(int2 p) { return p * 2 + int2(p.y & 1, p.x & 1); }

// ---- Ambient occlusion --------------------------------------------------------
// Scale free: no world-space radius. Samples spread log-uniformly between
// AO_R_MIN pixels and a tenth of the view's height, each judged at its own
// scale (its falloff radius AO_FALLOFF_SCALE times its distance), so the
// test is angular and the same constants hold at any scale of scene and
// camera; occluders far nearer the camera than they are beside the point fall
// outside every sample's falloff (no halos at silhouettes). A point fetch a
// sample from a level of the linear depth chosen by its distance (McGuire et
// al. 2012, "Scalable Ambient Obscurance"), each texel a real pixel's depth.

static const uint  AO_NUM_SAMPLES   = 16u;
static const int   AO_LOG_Q         = 3;      // a sample s pixels out reads level floor(log2 s) - AO_LOG_Q (1 at least)
static const int   AO_MAX_MIP       = 5;
static const float AO_FALLOFF_SCALE = 2.5;    // occluders steeper than ~66 degrees above a sample's scale fade out
static const float AO_BIAS          = 0.05;   // against self occlusion from depth quantisation
static const float AO_R_MIN         = 2.0;    // pixels
static const float AO_R_MAX         = 0.1;    // of the view's height
static const float GOLDEN_ANGLE     = 2.39996323;
static const int   BLUR_RADIUS      = 4;      // two periods of the 4 x 4 pattern
static const float BLUR_TOLERANCE   = 2.0;    // tangent-plane distance, pixel footprints

// Level `a.level` of the linear depth: 0 from md_render's depth, each other
// from the one before, a pixel kept per texel (the rotated grid)
[shader("compute")]
[numthreads(8, 8, 1)]
void ao_linear(uint3 tid : SV_DispatchThreadID) {
    Args  a = MD_ARGS;
    uint2 s = level_size(a, a.level);
    if (tid.x >= s.x || tid.y >= s.y) return;
    float z;
    if (a.level == 0u) {
        z = linear_depth(a, a.depth.Load(int3(int2(tid.xy), 0)));
    } else {
        z = lin_at(a, a.level - 1u, pick(int2(tid.xy)));
    }
    a.lin[level_offset(a, a.level) + tid.y * s.x + tid.x] = z;
}

static const uint BAYER[16] = { 0, 8, 2, 10, 12, 4, 14, 6, 3, 11, 1, 9, 15, 7, 13, 5 };

uint bit_reverse4(uint k) { return reversebits(k) >> 28u; }

// At half resolution: plane 0
[shader("compute")]
[numthreads(8, 8, 1)]
void ao_estimate(uint3 tid : SV_DispatchThreadID) {
    Args  a  = MD_ARGS;
    uint2 hs = level_size(a, 1u);
    if (tid.x >= hs.x || tid.y >= hs.y) return;
    int2  hp = int2(tid.xy);
    uint  o  = tid.y * hs.x + tid.x;
    float z  = lin_at(a, 1u, hp);
    if (z >= Z_MAX) { a.ao[o] = float2(1.0, z); return; }
    int2   fp   = pick(hp);
    float2 fpx  = float2(fp) + 0.5;
    float3 P    = view_at(a, fpx, z);
    float3 N    = normal_at(a, fp);
    float  pw   = px_world(a, z);
    float2 full = float2(float(a.width), float(a.height));
    float  r_max = max(AO_R_MAX * full.y, 2.0 * AO_R_MIN);

    // a 4 x 4 interleaved pattern: 16 (rotation, radial phase) pairs, Bayer
    // ordered; the blur integrates its periods
    uint   k   = BAYER[(hp.y & 3) * 4 + (hp.x & 3)];
    float  rk  = float(bit_reverse4(k)) / 16.0;
    float  rot = 6.2831853 * (float(k) + 0.5) / 16.0;
    float  jr  = frac(rk + 1.0 / 32.0);
    float2 cs  = float2(cos(rot), sin(rot));
    float  log_ratio = log2(r_max / AO_R_MIN);
    float  ao  = 0.0;
    for (uint i = 0u; i < AO_NUM_SAMPLES; ++i) {
        float  t  = (float(i) + frac(jr + float(i) * 0.618034)) / float(AO_NUM_SAMPLES);
        float  s  = AO_R_MIN * exp2(t * log_ratio);
        float  an = float(i) * GOLDEN_ANGLE;
        float2 d  = float2(cos(an), sin(an));
        d = float2(d.x * cs.x - d.y * cs.y, d.x * cs.y + d.y * cs.x);
        float2 spx = fpx + d * s;
        if (any(spx < float2(0.0)) || any(spx >= full)) continue;   // off screen: unoccluded
        int  m  = clamp(int(log2(s)) - AO_LOG_Q, 1, AO_MAX_MIP);
        int2 tx = int2(spx) >> m;
        uint2 ms = level_size(a, uint(m));
        tx = clamp(tx, int2(0), int2(ms) - 1);
        float sz = lin_at(a, uint(m), tx);
        if (sz >= Z_MAX) continue;
        // back down the picks to the pixel the depth is of, so that the point
        // is on the surface
        for (int l = m; l > 0; --l) tx = pick(tx);
        float3 v  = view_at(a, float2(tx) + 0.5, sz) - P;
        float  vv = dot(v, v);
        float  r  = AO_FALLOFF_SCALE * s * pw;
        float  fall = saturate(1.0 - vv / (r * r));
        float  vn = dot(v, N) * rsqrt(vv + 1e-12);
        ao += max(vn - AO_BIAS, 0.0) * fall;
    }
    ao *= 1.0 / (float(AO_NUM_SAMPLES) * (1.0 - AO_BIAS));
    a.ao[o] = float2(saturate(1.0 - a.pc.w * ao), z);
}

// The QUALITY estimator (viamd's ssao_bitmask.frag): scale-free horizon-based
// occlusion with a visibility bitmask, at half resolution, into plane 0 as
// ao_estimate. Each pixel walks AO_NUM_SLICES screen-space slices through
// itself (Jimenez et al. 2016); every depth sample along a slice occludes the
// arc between the directions to its front and to the same point pushed back
// by an assumed thickness (AO_THICKNESS times its lateral distance); the
// union of the arcs is a 32-sector bitmask per slice (Therrien et al. 2023),
// the sectors spaced to carry equal shares of the cosine-weighted
// hemisphere, so the visibility is 1 - occupied / 32. A thin occluder hides
// only the arc it covers, not everything below it. Out: the visibility to
// the power of the intensity x AO_POWER_SCALE (viamd's default intensity, 5,
// gives the plain cosine-weighted visibility).

static const uint  AO_NUM_SLICES  = 2u;
static const uint  AO_NUM_STEPS   = 8u;      // per side of a slice
static const float AO_THICKNESS   = 2.0;
static const float AO_POWER_SCALE = 0.2;
static const float HALF_PI        = 1.57079633;

// acos, Abramowitz & Stegun 4.4.45: |error| <= 6.7e-5
float fast_acos(float x) {
    float t = abs(x);
    float r = sqrt(1.0 - t) * (1.5707288 + t * (-0.2121144 + t * (0.0742610 + t * -0.0187293)));
    return x >= 0.0 ? r : PI - r;
}

// The cosine-weighted measure of the slice from the view direction (theta 0)
// to the in-plane direction (a, b) = (cos theta, sin theta), for a projected
// normal at angle gamma: sign(theta) (cos g + 2 theta sin g - cos(2 theta - g)) / 4
float slice_measure(float a, float b, float cos_g, float sin_g) {
    float theta    = (b < 0.0 ? -1.0 : 1.0) * fast_acos(a);
    float cos_2t_g = (a * a - b * b) * cos_g + 2.0 * a * b * sin_g;
    return (b < 0.0 ? -0.25 : 0.25) * (cos_g + 2.0 * theta * sin_g - cos_2t_g);
}

[shader("compute")]
[numthreads(8, 8, 1)]
void ao_estimate_quality(uint3 tid : SV_DispatchThreadID) {
    Args  a  = MD_ARGS;
    uint2 hs = level_size(a, 1u);
    if (tid.x >= hs.x || tid.y >= hs.y) return;
    int2  hp = int2(tid.xy);
    uint  o  = tid.y * hs.x + tid.x;
    float z  = lin_at(a, 1u, hp);
    if (z >= Z_MAX) { a.ao[o] = float2(1.0, z); return; }
    int2   fp   = pick(hp);
    float2 fpx  = float2(fp) + 0.5;
    float3 P    = view_at(a, fpx, z);
    float3 N    = normal_at(a, fp);
    bool   orth = a.misc.w > 0.5;
    float3 V    = orth ? float3(0.0, 0.0, 1.0) : -normalize(P);
    float  pw   = px_world(a, z);
    float2 full = float2(float(a.width), float(a.height));
    float  r_max = max(AO_R_MAX * full.y, 2.0 * AO_R_MIN);

    // the 4 x 4 interleaved (rotation, radial phase) pattern of ao_estimate
    uint   k    = BAYER[(hp.y & 3) * 4 + (hp.x & 3)];
    float  rk   = float(bit_reverse4(k)) / 16.0;
    float  phi0 = PI * (float(k) + 0.5) / (16.0 * float(AO_NUM_SLICES));
    float2 d    = float2(cos(phi0), sin(phi0));
    const float  SLICE_STEP = PI / float(AO_NUM_SLICES);
    const float2 slice_rot  = float2(cos(SLICE_STEP), sin(SLICE_STEP));

    // steps at r_min * ratio^(i + phase): an exp2 a side, then a multiply a step
    float log_ratio  = log2(r_max / AO_R_MIN);
    float step_ratio = exp2(log_ratio / float(AO_NUM_STEPS));
    float s_first0   = AO_R_MIN * exp2(rk * log_ratio / float(AO_NUM_STEPS));
    float s_first1   = AO_R_MIN * exp2(frac(rk + 0.5) * log_ratio / float(AO_NUM_STEPS));

    float vis_sum = 0.0, w_sum = 0.0;
    for (uint sl = 0u; sl < AO_NUM_SLICES; ++sl) {
        if (sl > 0u) d = float2(d.x * slice_rot.x - d.y * slice_rot.y, d.x * slice_rot.y + d.y * slice_rot.x);
        // the slice's frame: o in its plane, perpendicular to V, toward +d on
        // the screen (pixels run y down, view space y up)
        float3 dir3  = float3(d.x, -d.y, 0.0);
        float3 ov    = normalize(dir3 - V * dot(dir3, V));
        float3 axis  = cross(ov, V);
        float3 n_p   = N - axis * dot(N, axis);
        float  n_len = length(n_p);
        if (n_len < 1e-4) continue;
        float cos_g = clamp(dot(n_p, V) / n_len, -1.0, 1.0);
        float sgn_g = dot(n_p, ov) < 0.0 ? -1.0 : 1.0;
        float gamma = sgn_g * fast_acos(cos_g);
        float sin_g = sgn_g * sqrt(1.0 - cos_g * cos_g);
        float total = cos_g + gamma * sin_g;                        // measure(hi) - measure(lo)
        float m_lo  = -0.5 * (cos_g + (gamma - HALF_PI) * sin_g);   // measure(lo), lo = gamma - pi / 2
        float inv_total = 1.0 / total;
        // the arcs onto the sectors with a stochastic rounding: small arcs unbiased on average
        float bit_ofs = frac(rk * 7.0 + float(sl) * 0.618034);
        uint  bits = 0u;
        for (uint side = 0u; side < 2u; ++side) {
            float2 sd = side == 0u ? d : -d;
            float  s  = side == 0u ? s_first0 : s_first1;
            for (uint i = 0u; i < AO_NUM_STEPS; ++i, s *= step_ratio) {
                float2 spx = fpx + sd * s;
                if (any(spx < float2(0.0)) || any(spx >= full)) break;
                int   m  = clamp(int(firstbithigh(uint(s))) - AO_LOG_Q, 1, AO_MAX_MIP);
                int2  tx = clamp(int2(spx) >> m, int2(0), int2(level_size(a, uint(m))) - 1);
                float sz = lin_at(a, uint(m), tx);
                if (sz >= Z_MAX) continue;
                for (int l = m; l > 0; --l) tx = pick(tx);
                float3 S = view_at(a, float2(tx) + 0.5, sz);
                // its front, and the same point pushed back along its view ray
                float  thick = AO_THICKNESS * s * pw;
                float3 push  = orth ? float3(0.0, 0.0, -thick) : normalize(S) * thick;
                float3 Df = S - P;
                float  hf = dot(Df, n_p);
                float  hb = hf + dot(push, n_p);
                if (hf <= 0.0 && hb <= 0.0) continue;   // both behind the tangent plane
                float3 Db = Df + push;
                // their directions in the slice; outside the hemisphere, its edge on their side
                float2 ef = float2(dot(Df, V), dot(Df, ov));
                float2 eb = float2(dot(Db, V), dot(Db, ov));
                ef *= rsqrt(max(dot(ef, ef), 1e-30));
                eb *= rsqrt(max(dot(eb, eb), 1e-30));
                float uf = hf > 0.0 ? (slice_measure(ef.x, ef.y, cos_g, sin_g) - m_lo) * inv_total : (ef.y > 0.0 ? 1.0 : 0.0);
                float ub = hb > 0.0 ? (slice_measure(eb.x, eb.y, cos_g, sin_g) - m_lo) * inv_total : (eb.y > 0.0 ? 1.0 : 0.0);
                uint  b0 = uint(clamp(floor(min(uf, ub) * 32.0 + bit_ofs), 0.0, 32.0));
                uint  b1 = uint(clamp(floor(max(uf, ub) * 32.0 + bit_ofs), 0.0, 32.0));
                if (b1 > b0) {
                    uint count = b1 - b0;
                    bits |= (count >= 32u ? 0xFFFFFFFFu : ((1u << count) - 1u)) << b0;
                }
            }
        }
        // a slice's weight: the projected normal's length x its share of the hemisphere
        float w = n_len * total;
        vis_sum += w * (1.0 - float(countbits(bits)) / 32.0);
        w_sum   += w;
    }
    float vis = w_sum > 0.0 ? vis_sum / w_sum : 1.0;
    a.ao[o] = float2(pow(saturate(vis), a.pc.w * AO_POWER_SCALE), z);
}

// A separable bilateral blur at half resolution: a box with half-weight end
// taps over two periods of the pattern, stopped by the distance of a tap from
// the centre's tangent plane in pixel footprints (scale free)
[shader("compute")]
[numthreads(8, 8, 1)]
void ao_blur(uint3 tid : SV_DispatchThreadID) {
    Args  a  = MD_ARGS;
    uint2 hs = level_size(a, 1u);
    if (tid.x >= hs.x || tid.y >= hs.y) return;
    uint   n   = hs.x * hs.y;
    float2* sp = a.ao + a.src * n;
    float2* dp = a.ao + (1u - a.src) * n;
    int2   p   = int2(tid.xy);
    float2 c   = sp[tid.y * hs.x + tid.x];
    if (c.y >= Z_MAX) { dp[tid.y * hs.x + tid.x] = c; return; }
    int2   fp  = pick(p);
    float3 N   = normal_at(a, fp);
    float3 P   = view_at(a, float2(fp) + 0.5, c.y);
    float  inv_tol = 1.0 / (BLUR_TOLERANCE * px_world(a, c.y));
    int2   dir = a.dir == 0u ? int2(1, 0) : int2(0, 1);
    float  sum = c.x, wsum = 1.0;
    for (int i = -BLUR_RADIUS; i <= BLUR_RADIUS; ++i) {
        if (i == 0) continue;
        int2   q = clamp(p + dir * i, int2(0), int2(hs) - 1);
        float2 s = sp[uint(q.y) * hs.x + uint(q.x)];
        float3 S = view_at(a, float2(pick(q)) + 0.5, s.y);
        float  d = abs(dot(S - P, N)) * inv_tol;
        float  w = (abs(i) == BLUR_RADIUS ? 0.5 : 1.0) * saturate(1.0 - d);
        sum  += s.x * w;
        wsum += w;
    }
    dp[tid.y * hs.x + tid.x] = float2(sum / wsum, c.y);
}

// The occlusion at a full resolution pixel of linear depth z: a joint
// bilateral upsample of plane 0's 2 x 2 about it
float ao_at(Args a, int2 q, float z) {
    uint2  hs   = level_size(a, 1u);
    float2 hq   = (float2(q) + 0.5) * 0.5 - 0.5;
    float2 base = floor(hq);
    float2 f    = hq - base;
    float  fpw  = px_world(a, z);
    float  num = 0.0, den = 0.0;
    [ForceUnroll] for (int j = 0; j < 4; ++j) {
        int2   o  = int2(j & 1, j >> 1);
        int2   t  = clamp(int2(base) + o, int2(0), int2(hs) - 1);
        float2 s  = a.ao[uint(t.y) * hs.x + uint(t.x)];
        float  wb = (o.x == 1 ? f.x : 1.0 - f.x) * (o.y == 1 ? f.y : 1.0 - f.y);
        float  wz = 1.0 / (1.0 + abs(s.y - z) / (2.0 * fpw));
        float  w  = wb * wz * wz + 1e-5;
        num += s.x * w;
        den += w;
    }
    return num / den;
}

// ---- Shading ------------------------------------------------------------------

float3 srgb_to_linear(float3 c) {
    c = saturate(c);
    return lerp(c / 12.92, pow((c + 0.055) / 1.055, float3(2.4)), step(float3(0.04045), c));
}

float FresnelSchlick(float cosTheta, float F0) {
    return F0 + (1.0 - F0) * pow(saturate(1.0 - cosTheta), 5.0);
}

float DistributionGGX(float NdotH, float roughness) {
    float a     = roughness * roughness;
    float a2    = a * a;
    float denom = NdotH * NdotH * (a2 - 1.0) + 1.0;
    return a2 / (PI * denom * denom);
}

float GeometrySchlickGGX(float NdotV, float roughness) {
    float r = roughness + 1.0;
    float k = (r * r) / 8.0;
    return NdotV / (NdotV * (1.0 - k) + k);
}

float GeometrySmith(float NdotV, float NdotL, float roughness) {
    return GeometrySchlickGGX(NdotL, roughness) * GeometrySchlickGGX(NdotV, roughness);
}

// The directional albedo of the specular lobe, F0 * x + y (Karis 2014)
float2 EnvBRDFApprox(float NdotV, float roughness) {
    const float4 c0 = float4(-1.0, -0.0275, -0.572, 0.022);
    const float4 c1 = float4( 1.0,  0.0425,  1.04, -0.04);
    float4 r = roughness * c0 + c1;
    float  a004 = min(r.x * r.x, exp2(-9.28 * NdotV)) * r.x + r.y;
    return float2(-1.04, 1.04) * a004 + r.zw;
}

float3 environment_irradiance(Args a, float3 N) {
    return lerp(a.ground.xyz, a.sky.xyz, 0.5 + 0.5 * N.y);
}

float3 environment_radiance(Args a, float3 R, float roughness) {
    float w = 0.1 + roughness;
    return lerp(a.ground.xyz, a.sky.xyz, smoothstep(-w, w, R.y));
}

// The occlusion of the specular lobe from the ambient occlusion (Lagarde and
// de Rousiers 2014); roughness: alpha. 1 where ao is 1.
float specular_occlusion(float NdotV, float ao, float roughness) {
    return saturate(pow(NdotV + ao, exp2(-16.0 * roughness - 1.0)) - 1.0 + ao);
}

// A surface lit: what it reflects diffusely and specularly, and the part of
// the environment its specular lobe reflects (Fe). An opaque surface is their
// sum; a transparent one of coverage a (isosurface.frag's, isomesh.frag's
// shade_surface) reflects a x diffuse + specular, so a clear surface still
// shows its reflections, and lets 1 - opacity through,
// opacity = 1 - (1 - a)(1 - Fe).
struct Lit { float3 diffuse; float3 specular; float Fe; };

// compose_deferred.frag's shade, without shadows
Lit shade_pbr(Args a, float3 color, float3 V, float3 N, float ao) {
    float3 L = a.light.xyz;
    float3 H = normalize(L + V);
    float HdotV = saturate(dot(H, V));
    float NdotH = saturate(dot(N, H));
    float NdotL = saturate(dot(N, L));
    float NdotV = saturate(dot(N, V));
    float F0 = saturate(a.dir_radiance.w);
    float roughness = clamp(a.sky.w, 0.04, 1.0);
    float3 albedo = a.misc.z > 0.5 ? saturate(color) : srgb_to_linear(color);
    Lit r;
    r.diffuse = float3(0.0);
    r.specular = float3(0.0);
    if (NdotL > 0.0) {
        float  NDF = DistributionGGX(NdotH, roughness);
        float  G   = GeometrySmith(NdotV, NdotL, roughness);
        float  F   = FresnelSchlick(HdotV, F0);
        float3 E   = a.dir_radiance.xyz * a.light.w * NdotL;
        r.specular += float3(NDF * G * F / (4.0 * NdotV * NdotL + 0.0001)) * E;
        r.diffuse  += (1.0 - F) * albedo * ONE_OVER_PI * E;
    }
    float2 ab = EnvBRDFApprox(NdotV, roughness);
    r.Fe = F0 * ab.x + ab.y;
    r.diffuse  += (1.0 - r.Fe) * albedo * environment_irradiance(a, N) * ao;
    r.specular += r.Fe * environment_radiance(a, reflect(-V, N), roughness) * specular_occlusion(NdotV, ao, roughness * roughness);
    return r;
}

static const float BASIC_SHININESS = 48.0;

// compose_deferred.frag's shade_basic, without shadows; the occlusion
// darkens all of it. All diffuse: a transparent surface is its colour over
// its coverage.
Lit shade_basic(Args a, float3 color, float3 V, float3 N, float ao) {
    float3 L = a.light.xyz;
    float  w = 0.5 + 0.5 * dot(N, L);
    float  light = a.light.w;
    float  amb = a.ground.w;
    float3 c = saturate(color) * (amb + light * (1.0 - amb) * w * w);
    c /= max(1.0, max(c.r, max(c.g, c.b)));
    float3 H = normalize(L + V);
    float  s = light * a.misc.x * pow(max(dot(N, H), 0.0), BASIC_SHININESS);
    c = lerp(c, float3(1.0), saturate(s));
    Lit r;
    r.diffuse  = (a.misc.z > 0.5 ? c : srgb_to_linear(c)) * ao;
    r.specular = float3(0.0);
    r.Fe       = 0.0;
    return r;
}

Lit shade(Args a, float3 color, float3 V, float3 N, float ao) {
    return a.misc.y > 0.5 ? shade_basic(a, color, V, N, ao) : shade_pbr(a, color, V, N, ao);
}


// ---- Volumes ------------------------------------------------------------------
// A volume's record: VOL_STRIDE float4
//   0..3   view space to its texture space (columns)
//   4      clip_min, w: the transfer function's range begin
//   5      clip_max, w: 1 / its extent
//   6      its voxels per axis, w: 1 with a transfer function
//   7      x: isovalues, y: planes
//   8, 9   isovalues; 10..17 their colours; 18, 19 their optical densities
//   20..27 planes (texture space, n . x + d >= 0 kept)
//   28     the colour field's range: begin, 1 / extent; z: flags (1 colour
//          volume, 2 colour field, 4 normals); w: the isolines' width (pixels)
//   29     x: isolines, y: 1 with the transient one, z: its value
//   30, 31 the isolines' values
//   32..34 world to view (columns), for the normals
//   35     the colour field's voxels per axis

static const uint  VOL_STRIDE       = 36u;
static const uint  MAX_VOL          = 4u;
static const uint  MAX_ISO          = 8u;
static const uint  MAX_STEPS        = 4096u;
static const float REF_SAMPLING     = 150.0;   // a transfer function's opacity is of a slab 1/150 of the box
static const float STEPS_PER_VOXEL  = 2.0;
static const float T_MIN            = 0.005;   // below this transmittance the ray is opaque
static const float T_FAR            = 1.0e30;

float4 vrec(Args a, uint v, uint i) { return a.vols[v * VOL_STRIDE + i]; }

float vol_field(Args a, uint v, float3 p) {
    switch (v) {
    case 0u: return a.vf0.SampleLevel(a.trilinear, p, 0.0);
    case 1u: return a.vf1.SampleLevel(a.trilinear, p, 0.0);
    case 2u: return a.vf2.SampleLevel(a.trilinear, p, 0.0);
    default: return a.vf3.SampleLevel(a.trilinear, p, 0.0);
    }
}

float4 vol_transfer(Args a, uint v, float t) {
    float2 uv = float2(t, 0.5);
    switch (v) {
    case 0u: return a.vt0.SampleLevel(a.trilinear, uv, 0.0);
    case 1u: return a.vt1.SampleLevel(a.trilinear, uv, 0.0);
    case 2u: return a.vt2.SampleLevel(a.trilinear, uv, 0.0);
    default: return a.vt3.SampleLevel(a.trilinear, uv, 0.0);
    }
}

float4 vol_color(Args a, uint v, float3 p) {
    switch (v) {
    case 0u: return a.vc0.SampleLevel(a.trilinear, p, 0.0);
    case 1u: return a.vc1.SampleLevel(a.trilinear, p, 0.0);
    case 2u: return a.vc2.SampleLevel(a.trilinear, p, 0.0);
    default: return a.vc3.SampleLevel(a.trilinear, p, 0.0);
    }
}

float vol_cfield(Args a, uint v, float3 p) {
    switch (v) {
    case 0u: return a.vg0.SampleLevel(a.trilinear, p, 0.0);
    case 1u: return a.vg1.SampleLevel(a.trilinear, p, 0.0);
    case 2u: return a.vg2.SampleLevel(a.trilinear, p, 0.0);
    default: return a.vg3.SampleLevel(a.trilinear, p, 0.0);
    }
}

float4 vol_cmap(Args a, uint v, float t) {
    float2 uv = float2(t, 0.5);
    switch (v) {
    case 0u: return a.vm0.SampleLevel(a.trilinear, uv, 0.0);
    case 1u: return a.vm1.SampleLevel(a.trilinear, uv, 0.0);
    case 2u: return a.vm2.SampleLevel(a.trilinear, uv, 0.0);
    default: return a.vm3.SampleLevel(a.trilinear, uv, 0.0);
    }
}

float4 vol_normals(Args a, uint v, float3 p) {
    switch (v) {
    case 0u: return a.vn0.SampleLevel(a.trilinear, p, 0.0);
    case 1u: return a.vn1.SampleLevel(a.trilinear, p, 0.0);
    case 2u: return a.vn2.SampleLevel(a.trilinear, p, 0.0);
    default: return a.vn3.SampleLevel(a.trilinear, p, 0.0);
    }
}

float iso_value(Args a, uint v, uint i)   { float4 q = vrec(a, v, 8u + i / 4u);  return q[i % 4u]; }
float iso_density(Args a, uint v, uint i) { float4 q = vrec(a, v, 18u + i / 4u); return q[i % 4u]; }

bool is_inside(float f, float iv) { return iv >= 0.0 ? f >= iv : f <= iv; }

uint membership(Args a, uint v, float f) {
    uint m = 0u, n = uint(vrec(a, v, 7u).x);
    for (uint i = 0u; i < n; ++i) if (is_inside(f, iso_value(a, v, i))) m |= 1u << i;
    return m;
}

float tau_of(Args a, uint v, uint m) {
    float t = 0.0;
    uint  n = uint(vrec(a, v, 7u).x);
    for (uint i = 0u; i < n; ++i) if ((m & (1u << i)) != 0u) t += max(iso_density(a, v, i), 0.0);
    return t;
}

// Where between ta and tb (fa, fb the field there, on either side of iv)
// the field along o + t d takes the value iv: false position and bisection
// in turn, on the trilinear field
float refine(Args a, uint v, float3 o, float3 d, float ta, float tb, float fa, float fb, float iv) {
    float ga = fa - iv, gb = fb - iv;
    for (int k = 0; k < 6; ++k) {
        float t = (k & 1) == 0 && ga != gb ? clamp(ta + (tb - ta) * ga / (ga - gb), ta, tb) : 0.5 * (ta + tb);
        float g = vol_field(a, v, o + d * t) - iv;
        if ((g < 0.0) == (ga < 0.0)) { ta = t; ga = g; } else { tb = t; gb = g; }
    }
    return ga != gb ? ta + (tb - ta) * ga / (ga - gb) : 0.5 * (ta + tb);
}

// A texture-space gradient as a view-space vector: the transpose of view to
// texture's linear part
float3 grad_to_view(Args a, uint v, float3 g) {
    return float3(dot(vrec(a, v, 0u).xyz, g), dot(vrec(a, v, 1u).xyz, g), dot(vrec(a, v, 2u).xyz, g));
}

// The field's gradient at p (texture space) as a view-space direction: the
// normals volume's where it has one not too short, else central differences
// a voxel apart
float3 vol_normal(Args a, uint v, float3 p) {
    if ((uint(vrec(a, v, 28u).z) & 4u) != 0u) {
        float3 n = vol_normals(a, v, p).xyz * 2.0 - 1.0;
        if (dot(n, n) > 0.25) return vrec(a, v, 32u).xyz * n.x + vrec(a, v, 33u).xyz * n.y + vrec(a, v, 34u).xyz * n.z;
    }
    float3 dim = vrec(a, v, 6u).xyz;
    float3 h   = 1.0 / dim;
    float3 g   = float3(vol_field(a, v, p + float3(h.x, 0, 0)) - vol_field(a, v, p - float3(h.x, 0, 0)),
                        vol_field(a, v, p + float3(0, h.y, 0)) - vol_field(a, v, p - float3(0, h.y, 0)),
                        vol_field(a, v, p + float3(0, 0, h.z)) - vol_field(a, v, p - float3(0, 0, h.z))) * dim;
    return grad_to_view(a, v, g);
}

// ---- Isolines (viamd's isolines.h, isosurface.frag) ------------------------------
// A line of value u is |f - u| / g pixels away, g the field's change per
// pixel along the surface: its gradient along where the rays of the next
// pixels across and down meet the surface's tangent plane (ray differentials)

float isoline_coverage(float d, float w) { return saturate(0.5 * w + 0.5 - d); }

float4 apply_isolines(Args a, uint v, float4 c, float f, float g) {
    const float3 DARK  = float3(0.02);
    const float3 LIGHT = float3(0.95);
    const float  RIM   = 2.0;           // pixels, both sides together
    float4 r  = vrec(a, v, 28u), h = vrec(a, v, 29u);
    float  w  = r.w, ig = 1.0 / max(g, 1.0e-30);
    float  core = 0.0, rim = 0.0;
    uint   n = uint(h.x);
    for (uint i = 0u; i < n; ++i) {
        float4 q = vrec(a, v, 30u + i / 4u);
        float  d = abs(f - q[i % 4u]) * ig;
        core = max(core, isoline_coverage(d, w));
        rim  = max(rim,  isoline_coverage(d, w + RIM));
    }
    c.rgb = lerp(c.rgb, LIGHT, rim);
    c.rgb = lerp(c.rgb, DARK, core);
    float cover = rim;
    if (h.y > 0.5) {
        float d  = abs(f - h.z) * ig;
        float hc = isoline_coverage(d, w), hr = isoline_coverage(d, w + RIM);
        c.rgb = lerp(c.rgb, DARK, hr);
        c.rgb = lerp(c.rgb, LIGHT, hc);
        cover = max(cover, hr);
    }
    c.a = lerp(c.a, 1.0, cover);
    return c;
}

// Where the ray o + s d meets the plane through P of normal N (at a grazing
// angle far off: the field changes a lot from pixel to pixel there)
float3 plane_hit(float3 o, float3 d, float3 P, float3 N) {
    float dn = dot(N, d);
    dn = dn < 0.0 ? min(dn, -1.0e-6) : max(dn, 1.0e-6);
    return o + d * (dot(N, P - o) / dn);
}

// The range of t (o + t d, texture space) inside the volume's clip box and
// on the kept side of its planes, within [lo, hi]
float2 vol_range(Args a, uint v, float3 o, float3 d, float lo, float hi) {
    float3 bmin = vrec(a, v, 4u).xyz, bmax = vrec(a, v, 5u).xyz;
    float3 inv  = 1.0 / select(abs(d) < float3(1e-20), float3(1e-20), d);
    float3 t0 = (bmin - o) * inv, t1 = (bmax - o) * inv;
    float3 tn = min(t0, t1), tf = max(t0, t1);
    lo = max(lo, max(tn.x, max(tn.y, tn.z)));
    hi = min(hi, min(tf.x, min(tf.y, tf.z)));
    uint np = uint(vrec(a, v, 7u).y);
    for (uint i = 0u; i < np; ++i) {
        float4 pl = vrec(a, v, 20u + i);
        float  pa = dot(pl.xyz, o) + pl.w, pb = dot(pl.xyz, d);
        if (abs(pb) < 1e-20) { if (pa < 0.0) hi = lo - 1.0; continue; }
        float t = -pa / pb;
        if (pb > 0.0) lo = max(lo, t); else hi = min(hi, t);
    }
    return float2(lo, hi);
}

uint hash_u(uint x) {
    x ^= x >> 16; x *= 0x7feb352du; x ^= x >> 15; x *= 0x846ca68bu; x ^= x >> 16;
    return x;
}

// ---- The compose --------------------------------------------------------------

float3 unproject(Args a, float2 ndc, float d) {
    float4 h = a.ip0 * ndc.x + a.ip1 * ndc.y + a.ip2 * d + a.ip3;
    return h.xyz / h.w;
}

// Isosurface i of volume v's colour at pt (texture space), P (view space),
// of view-space normal N, seen through the pixel at ndc
float4 surface_color(Args a, uint v, uint i, float3 pt, float3 P, float3 N, float2 ndc) {
    float4 c = vrec(a, v, 10u + i);
    float4 r = vrec(a, v, 28u);
    uint   fl = uint(r.z);
    if ((fl & 1u) != 0u) return c * vol_color(a, v, pt);
    if ((fl & 2u) == 0u) return c;
    float f = vol_cfield(a, v, pt);
    c.rgb *= vol_cmap(a, v, saturate((f - r.x) * r.y)).rgb;
    float4 h = vrec(a, v, 29u);
    if (h.x < 0.5 && h.y < 0.5) return c;
    // its change per pixel: the gradient (central differences a voxel of the
    // colour field apart) along where the next pixels' rays meet the tangent plane
    float3 dim = vrec(a, v, 35u).xyz;
    float3 e   = 1.0 / dim;
    float3 g   = float3(vol_cfield(a, v, pt + float3(e.x, 0, 0)) - vol_cfield(a, v, pt - float3(e.x, 0, 0)),
                        vol_cfield(a, v, pt + float3(0, e.y, 0)) - vol_cfield(a, v, pt - float3(0, e.y, 0)),
                        vol_cfield(a, v, pt + float3(0, 0, e.z)) - vol_cfield(a, v, pt - float3(0, 0, e.z))) * dim * 0.5;
    float3 G   = grad_to_view(a, v, g);
    float2 nx  = ndc + float2(2.0 / float(a.width), 0.0);
    float2 ny  = ndc - float2(0.0, 2.0 / float(a.height));
    float3 ox  = unproject(a, nx, 1.0), oy = unproject(a, ny, 1.0);
    float3 dx  = a.misc.w > 0.5 ? float3(0.0, 0.0, -1.0) : normalize(ox);
    float3 dy  = a.misc.w > 0.5 ? float3(0.0, 0.0, -1.0) : normalize(oy);
    float  gp  = length(float2(dot(G, plane_hit(ox, dx, P, N) - P), dot(G, plane_hit(oy, dy, P, N) - P)));
    return apply_isolines(a, v, c, f, gp);
}

// ---- Slices -------------------------------------------------------------------
// A slice's record: SLICE_STRIDE float4
//   0..3   view space to its rectangle's (the unit square at z 0)
//   4..7   view space to its field's texture space
//   8      the transfer function's range begin, 1 / its extent, the opacity, the edge's half width (pixels)
//   9      the edge's colour (linear), alpha
//   10..13 its rectangle's space to view space

static const uint SLICE_STRIDE = 14u;
static const uint MAX_SLICES   = 8u;

float4 srec(Args a, uint i, uint k) { return a.slices[i * SLICE_STRIDE + k]; }

float slice_field(Args a, uint i, float3 p) {
    switch (i) {
    case 0u: return a.sf0.SampleLevel(a.trilinear, p, 0.0);
    case 1u: return a.sf1.SampleLevel(a.trilinear, p, 0.0);
    case 2u: return a.sf2.SampleLevel(a.trilinear, p, 0.0);
    case 3u: return a.sf3.SampleLevel(a.trilinear, p, 0.0);
    case 4u: return a.sf4.SampleLevel(a.trilinear, p, 0.0);
    case 5u: return a.sf5.SampleLevel(a.trilinear, p, 0.0);
    case 6u: return a.sf6.SampleLevel(a.trilinear, p, 0.0);
    default: return a.sf7.SampleLevel(a.trilinear, p, 0.0);
    }
}

float3 slice_transfer(Args a, uint i, float t) {
    float2 uv = float2(t, 0.5);
    switch (i) {
    case 0u: return a.st0.SampleLevel(a.trilinear, uv, 0.0).rgb;
    case 1u: return a.st1.SampleLevel(a.trilinear, uv, 0.0).rgb;
    case 2u: return a.st2.SampleLevel(a.trilinear, uv, 0.0).rgb;
    case 3u: return a.st3.SampleLevel(a.trilinear, uv, 0.0).rgb;
    case 4u: return a.st4.SampleLevel(a.trilinear, uv, 0.0).rgb;
    case 5u: return a.st5.SampleLevel(a.trilinear, uv, 0.0).rgb;
    case 6u: return a.st6.SampleLevel(a.trilinear, uv, 0.0).rgb;
    default: return a.st7.SampleLevel(a.trilinear, uv, 0.0).rgb;
    }
}

float3 mulp4(float4 c0, float4 c1, float4 c2, float4 c3, float3 p) { return (c0 * p.x + c1 * p.y + c2 * p.z + c3).xyz; }
float3 mulv4(float4 c0, float4 c1, float4 c2, float3 v) { return (c0 * v.x + c1 * v.y + c2 * v.z).xyz; }

// Where the ray o + t d (view space) meets slice i's plane: t, and the point
// on its rectangle in x, y; t < 0 when it does not
float3 slice_hit(Args a, uint i, float3 o, float3 d) {
    float3 ro = mulp4(srec(a, i, 0u), srec(a, i, 1u), srec(a, i, 2u), srec(a, i, 3u), o);
    float3 rd = mulv4(srec(a, i, 0u), srec(a, i, 1u), srec(a, i, 2u), d);
    if (abs(rd.z) < 1e-20) return float3(0.0, 0.0, -1.0);
    float t = -ro.z / rd.z;
    float2 q = (ro + rd * t).xy;
    return float3(q, t);
}

// The pixels from the pixel at ndc to the point b of slice i's rectangle, projected; far behind the camera
float slice_px(Args a, uint i, float2 b, float2 ndc) {
    float3 v = mulp4(srec(a, i, 10u), srec(a, i, 11u), srec(a, i, 12u), srec(a, i, 13u), float3(b, 0.0));
    float  w = a.pc.x * v.z + a.pc.y;
    if (w <= 0.0) return 1e9;
    float2 n = float2(a.pa.x * v.x + a.pa.z * v.z + a.pb.x, a.pa.y * v.y + a.pa.w * v.z + a.pb.y) / w;
    return length((n - ndc) * float2(float(a.width), float(a.height)) * 0.5);
}

// A film shaded onto the ray: what it reflects, and what it lets through
void film(Args a, inout float3 C, inout float T, float3 color, float alpha, float3 V, float3 N) {
    Lit l = shade(a, color, V, N, 1.0);
    C += T * (alpha * l.diffuse + l.specular);
    T *= (1.0 - alpha) * (1.0 - l.Fe);
}

// Slice i onto the ray at t, d pixels in from its rectangle's border (out:
// negative): its edge in front of it, the line over its halo, each feathered
// over a pixel; the slice within its rectangle, its border feathered too
void slice(Args a, inout float3 C, inout float T, uint i, float3 P, float d) {
    float4 r = srec(a, i, 8u), e = srec(a, i, 9u);
    float  w = r.w;
    if (e.a > 0.0) {
        float core = saturate(w + 0.5 - abs(d)) * e.a;
        float halo = saturate(w + 1.5 - abs(d)) * e.a * 0.45;
        C += T * core * e.rgb;
        T *= 1.0 - (core + halo * (1.0 - core));
    }
    float al = r.z * saturate(d + 0.5);
    if (al > 0.0) {
        float3 pf = mulp4(srec(a, i, 4u), srec(a, i, 5u), srec(a, i, 6u), srec(a, i, 7u), P);
        float  u  = saturate((slice_field(a, i, pf) - r.x) * r.y);
        C += T * al * srgb_to_linear(slice_transfer(a, i, u));
        T *= 1.0 - al;
    }
}

// Layer k of the pixel p onto the ray
void layer(Args a, inout float3 C, inout float T, uint k, uint p, float3 V) {
    uint   c8 = a.layers[k * a.stride + p];
    float4 c  = float4(float(c8 & 0xFFu), float((c8 >> 8) & 0xFFu), float((c8 >> 16) & 0xFFu), float(c8 >> 24)) / 255.0;
    uint   nb = a.layers[a.count * a.stride + k * a.stride + p];
    float3 N  = oct_decode(float2(half_to_float(nb & 0xFFFFu), half_to_float(nb >> 16)));
    film(a, C, T, c.rgb, c.a, V, N);
}

// The surface kd (a layer, or 0x100 | a slice) onto the ray at P
void surface(Args a, inout float3 C, inout float T, uint kd, float d, uint p, float3 P, float3 V) {
    if ((kd & 0x100u) != 0u) slice(a, C, T, kd & 0xFFu, P, d);
    else                     layer(a, C, T, kd, p, V);
}

// The events along the ray in a step, in order: (t, v | i << 8), i an
// isosurface crossed, 254 the volume v begins, 255 it ends. Parallel arrays
// of the kernel's own, filled in place (no arrays through parameters: Slang
// copies them in and out).
static const uint MAX_HITS = 12u;
static const uint EV_BEGIN = 254u;
static const uint EV_END   = 255u;

[shader("compute")]
[numthreads(8, 8, 1)]
void compose(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    if (tid.x >= a.width || tid.y >= a.height) return;
    int2   q   = int2(tid.xy);
    uint   p   = tid.y * a.width + tid.x;
    float2 ndc = ndc_of(a, float2(q) + 0.5);

    // the ray: from the near plane, toward the eye V
    float3 O = unproject(a, ndc, 1.0);
    float3 D = a.misc.w > 0.5 ? float3(0.0, 0.0, -1.0) : normalize(O);
    float3 V = -D;

    // the opaque surface
    float3 S   = float3(0.0);
    float4 col = a.color.Load(int3(q, 0));
    float  cov = col.a > 0.0 ? 1.0 : 0.0;
    float  dop = a.depth.Load(int3(q, 0));
    float  top = T_FAR;
    if (cov > 0.0) {
        float3 N  = normal_at(a, q);
        float  ao = a.use_ao != 0u ? ao_at(a, q, a.lin[tid.y * a.width + tid.x]) : 1.0;
        Lit l = shade(a, col.rgb, V, N, ao);
        S   = l.diffuse + l.specular;
        top = dot(unproject(a, ndc, dop) - O, D);
    }

    // the layers, front to back, and where they are along the ray
    // and the slices it meets in front of the opaque surface: the surfaces
    // along the ray, in order: (t, k) for layer k, (t, 0x100 | i) slice i;
    // a slice's distance in from its border (pixels) beside it
    uint  nl = 0u;
    float tl[16];
    uint  kd[16];
    float sd[16];
    for (uint k = 0u; k < a.count && k < 8u; ++k) {
        float d = asfloat(a.layers[2u * a.count * a.stride + k * a.stride + p]);
        if (d == 0.0) break;
        tl[k] = dot(unproject(a, ndc, d) - O, D);
        kd[k] = k;
        sd[k] = 0.0;
        nl = k + 1u;
    }
    uint ns = min(a.slice_count, MAX_SLICES);
    for (uint i = 0u; i < ns; ++i) {
        float3 h = slice_hit(a, i, O, D);
        if (h.z < 0.0 || h.z > top) continue;
        // the pixels in from the rectangle's border: to its nearest point projected, which holds where the plane
        // runs off to its horizon, as an estimate from the rays of the next pixels does not (the whole of the far
        // plane would be within a pixel of the border). Within it, the nearest of the points across on its sides.
        float d;
        if (all(h.xy >= 0.0) && all(h.xy <= 1.0)) {
            d = min(min(slice_px(a, i, float2(0.0, h.y), ndc), slice_px(a, i, float2(1.0, h.y), ndc)),
                    min(slice_px(a, i, float2(h.x, 0.0), ndc), slice_px(a, i, float2(h.x, 1.0), ndc)));
        } else {
            d = -slice_px(a, i, saturate(h.xy), ndc);
        }
        float4 r = srec(a, i, 8u), e = srec(a, i, 9u);
        if (!((r.z > 0.0 && d > -0.5) || (e.a > 0.0 && abs(d) < r.w + 1.5))) continue;
        uint j = nl++;
        while (j > 0u && tl[j - 1u] > h.z) { tl[j] = tl[j - 1u]; kd[j] = kd[j - 1u]; sd[j] = sd[j - 1u]; --j; }
        tl[j] = h.z; kd[j] = 0x100u | i; sd[j] = d;
    }

    float3 C = float3(0.0);
    float  T = 1.0;
    uint   kl = 0u;     // the next layer

    // the volumes: each ray in its texture space, the part of it inside
    uint   nv = min(a.vol_count, MAX_VOL);
    float3 vo[MAX_VOL], vd[MAX_VOL];
    float2 vr[MAX_VOL];
    float  t_beg = T_FAR, t_end = -T_FAR, rate = 0.0;
    for (uint v = 0u; v < nv; ++v) {
        float4 c0 = vrec(a, v, 0u), c1 = vrec(a, v, 1u), c2 = vrec(a, v, 2u), c3 = vrec(a, v, 3u);
        vo[v] = (c0 * O.x + c1 * O.y + c2 * O.z + c3).xyz;
        vd[v] = (c0 * D.x + c1 * D.y + c2 * D.z).xyz;
        vr[v] = vol_range(a, v, vo[v], vd[v], 0.0, top);
        if (vr[v].x < vr[v].y) {
            t_beg = min(t_beg, vr[v].x);
            t_end = max(t_end, vr[v].y);
            rate  = max(rate, length(vd[v] * vrec(a, v, 6u).xyz));     // voxels per unit of t
        }
    }
    if (t_beg < t_end && rate > 0.0) {
        uint   n   = min(uint(ceil((t_end - t_beg) * rate * STEPS_PER_VOXEL)), MAX_STEPS);
        n = max(n, 1u);
        float  dt  = (t_end - t_beg) / float(n);
        float  jit = float(hash_u(p * 0x9E3779B9u ^ hash_u(a.frame)) >> 8) / 16777216.0;
        float  fprev[MAX_VOL];
        uint   inside[MAX_VOL];
        uint   state = 0u;          // bit v: volume v entered; bit 8 + v: left
        float  tau = 0.0, t_abs = t_beg;
        [ForceUnroll] for (uint v = 0u; v < MAX_VOL; ++v) { fprev[v] = 0.0; inside[v] = 0u; }
        for (uint s = 0u; s < n && T >= T_MIN; ++s) {
            float ta = t_beg + dt * float(s), tb = ta + dt, ts = ta + dt * jit;
            // the crossings and the volumes' ends in (ta, tb], in order along the ray
            float ht[MAX_HITS];
            uint  he[MAX_HITS];
            uint  nh = 0u;
            for (uint v = 0u; v < nv; ++v) {
                if (vr[v].x >= vr[v].y || vr[v].y < ta || vr[v].x > tb || (state & (0x100u << v)) != 0u) continue;
                float ua = max(ta, vr[v].x), ub = min(tb, vr[v].y);
                if ((state & (1u << v)) == 0u) {     // it begins here: its densities from there on
                    fprev[v]  = vol_field(a, v, vo[v] + vd[v] * ua);
                    inside[v] = membership(a, v, fprev[v]);
                    state    |= 1u << v;
                    { float et = ua; uint ee = v | (EV_BEGIN << 8); if (nh < MAX_HITS) { uint j = nh++; while (j > 0u && ht[j - 1u] > et) { ht[j] = ht[j - 1u]; he[j] = he[j - 1u]; --j; } ht[j] = et; he[j] = ee; } }
                }
                float fb = vol_field(a, v, vo[v] + vd[v] * ub);
                uint  ni = uint(vrec(a, v, 7u).x);
                for (uint i = 0u; i < ni; ++i) {
                    float iv = iso_value(a, v, i);
                    if (is_inside(fprev[v], iv) == is_inside(fb, iv)) continue;
                    { float et = refine(a, v, vo[v], vd[v], ua, ub, fprev[v], fb, iv); uint ee = v | (i << 8); if (nh < MAX_HITS) { uint j = nh++; while (j > 0u && ht[j - 1u] > et) { ht[j] = ht[j - 1u]; he[j] = he[j - 1u]; --j; } ht[j] = et; he[j] = ee; } }
                }
                fprev[v] = fb;
                if (ub >= vr[v].y) { float et = vr[v].y; uint ee = v | (EV_END << 8); if (nh < MAX_HITS) { uint j = nh++; while (j > 0u && ht[j - 1u] > et) { ht[j] = ht[j - 1u]; he[j] = he[j - 1u]; --j; } ht[j] = et; he[j] = ee; } }    // it ends here: no more of its density
            }
            // the events in order: crossings, layers, the emission's sample
            uint ih = 0u;
            bool emitted = false;
            for (uint guard = 0u; guard < 32u; ++guard) {
                float th = ih < nh ? ht[ih] : T_FAR;
                float tk = kl < nl && tl[kl] <= tb ? tl[kl] : T_FAR;
                float te = emitted ? T_FAR : ts;
                float t  = min(th, min(tk, te));
                if (t >= T_FAR) break;
                if (tau > 0.0) T *= exp(-tau * max(t - t_abs, 0.0));
                t_abs = max(t_abs, t);
                if (t == tk) {
                    surface(a, C, T, kd[kl], sd[kl], p, O + D * tk, V);
                    ++kl;
                } else if (t == th) {
                    uint ev = he[ih] & 0xFFu, ei = he[ih] >> 8;
                    ++ih;
                    if (ei == EV_BEGIN) {
                        tau += tau_of(a, ev, inside[ev]);
                    } else if (ei == EV_END) {
                        tau -= tau_of(a, ev, inside[ev]);
                        inside[ev] = 0u;
                        state |= 0x100u << ev;
                    } else {
                        float3 pt = vo[ev] + vd[ev] * t;
                        float3 N  = vol_normal(a, ev, pt);
                        float  ln = length(N);
                        N = ln > 1e-20 ? N / ln : V;
                        if (dot(N, V) < 0.0) N = -N;
                        float4 c = surface_color(a, ev, ei, pt, O + D * t, N, ndc);
                        film(a, C, T, c.rgb, saturate(c.a), V, N);
                        tau -= tau_of(a, ev, inside[ev]);
                        inside[ev] ^= 1u << ei;
                        tau += tau_of(a, ev, inside[ev]);
                    }
                } else {
                    emitted = true;
                    // each volume there through its transfer function, the
                    // opacity for the step's length in texture space
                    for (uint v = 0u; v < nv; ++v) {
                        if (vrec(a, v, 6u).w < 0.5 || ts < vr[v].x || ts > vr[v].y) continue;
                        float  f  = vol_field(a, v, vo[v] + vd[v] * ts);
                        float  u  = saturate((f - vrec(a, v, 4u).w) * vrec(a, v, 5u).w);
                        float4 tf = vol_transfer(a, v, u);
                        float  w  = dt * length(vd[v]) * REF_SAMPLING;
                        float  al = 1.0 - pow(max(1.0 - tf.a, 1e-6), w);
                        C += T * al * srgb_to_linear(tf.rgb);
                        T *= 1.0 - al;
                    }
                }
            }
        }
        if (tau > 0.0 && t_end > t_abs) T *= exp(-tau * (t_end - t_abs));
    }
    // the layers behind the volumes
    for (; kl < nl; ++kl) surface(a, C, T, kd[kl], sd[kl], p, O + D * tl[kl], V);
    a.target[q] = float4(C + T * cov * S, 1.0 - T * (1.0 - cov));
}
