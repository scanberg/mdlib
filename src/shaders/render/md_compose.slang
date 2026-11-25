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

// compose_deferred.frag's shade, without shadows
float3 shade_pbr(Args a, float3 color, float3 V, float3 N, float ao) {
    float3 L = a.light.xyz;
    float3 H = normalize(L + V);
    float HdotV = saturate(dot(H, V));
    float NdotH = saturate(dot(N, H));
    float NdotL = saturate(dot(N, L));
    float NdotV = saturate(dot(N, V));
    float F0 = saturate(a.dir_radiance.w);
    float roughness = clamp(a.sky.w, 0.04, 1.0);
    float3 albedo = a.misc.z > 0.5 ? saturate(color) : srgb_to_linear(color);
    float3 Lo = float3(0.0);
    if (NdotL > 0.0) {
        float NDF = DistributionGGX(NdotH, roughness);
        float G   = GeometrySmith(NdotV, NdotL, roughness);
        float F   = FresnelSchlick(HdotV, F0);
        float3 specular = float3(NDF * G * F / (4.0 * NdotV * NdotL + 0.0001));
        Lo += ((1.0 - F) * albedo * ONE_OVER_PI + specular) * (a.dir_radiance.xyz * a.light.w) * NdotL;
    }
    float2 ab = EnvBRDFApprox(NdotV, roughness);
    float  Fe = F0 * ab.x + ab.y;
    float3 R  = reflect(-V, N);
    Lo += (1.0 - Fe) * albedo * environment_irradiance(a, N) * ao +
          Fe * environment_radiance(a, R, roughness) * specular_occlusion(NdotV, ao, roughness * roughness);
    return Lo;
}

static const float BASIC_SHININESS = 48.0;

// compose_deferred.frag's shade_basic, without shadows; the occlusion
// darkens all of it
float3 shade_basic(Args a, float3 color, float3 V, float3 N, float ao) {
    float3 L = a.light.xyz;
    float  w = 0.5 + 0.5 * dot(N, L);
    float  light = a.light.w;
    float  amb = a.ground.w;
    float3 c = saturate(color) * (amb + light * (1.0 - amb) * w * w);
    c /= max(1.0, max(c.r, max(c.g, c.b)));
    float3 H = normalize(L + V);
    float  s = light * a.misc.x * pow(max(dot(N, H), 0.0), BASIC_SHININESS);
    c = lerp(c, float3(1.0), saturate(s));
    return (a.misc.z > 0.5 ? c : srgb_to_linear(c)) * ao;
}

float3 shade(Args a, float3 color, float3 V, float3 N, float ao) {
    return a.misc.y > 0.5 ? shade_basic(a, color, V, N, ao) : shade_pbr(a, color, V, N, ao);
}

// Toward the eye from the point at a pixel's centre and depth d (reverse Z)
float3 toward_eye(Args a, float2 ndc, float d) {
    if (a.misc.w > 0.5) return float3(0.0, 0.0, 1.0);
    float4 h = a.ip0 * ndc.x + a.ip1 * ndc.y + a.ip2 * d + a.ip3;
    return -normalize(h.xyz / h.w);
}

// ---- The compose --------------------------------------------------------------

[shader("compute")]
[numthreads(8, 8, 1)]
void compose(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    if (tid.x >= a.width || tid.y >= a.height) return;
    int2   q   = int2(tid.xy);
    uint   p   = tid.y * a.width + tid.x;
    float2 ndc = ndc_of(a, float2(q) + 0.5);

    // the opaque surface
    float3 S   = float3(0.0);
    float4 col = a.color.Load(int3(q, 0));
    float  cov = col.a > 0.0 ? 1.0 : 0.0;
    if (cov > 0.0) {
        float  d  = a.depth.Load(int3(q, 0));
        float3 N  = normal_at(a, q);
        float  ao = a.use_ao != 0u ? ao_at(a, q, a.lin[tid.y * a.width + tid.x]) : 1.0;
        S = shade(a, col.rgb, toward_eye(a, ndc, d), N, ao);
    }

    // the layers over it, front to back
    float3 C = float3(0.0);
    float  T = 1.0;
    for (uint k = 0u; k < a.count; ++k) {
        float d = asfloat(a.layers[2u * a.count * a.stride + k * a.stride + p]);
        if (d == 0.0) break;
        uint   c8 = a.layers[k * a.stride + p];
        float4 c  = float4(float(c8 & 0xFFu), float((c8 >> 8) & 0xFFu), float((c8 >> 16) & 0xFFu), float(c8 >> 24)) / 255.0;
        uint   nb = a.layers[a.count * a.stride + k * a.stride + p];
        float3 N  = oct_decode(float2(half_to_float(nb & 0xFFFFu), half_to_float(nb >> 16)));
        float3 L  = shade(a, c.rgb, toward_eye(a, ndc, d), N, 1.0);
        C += T * c.a * L;
        T *= 1.0 - c.a;
    }
    a.target[q] = float4(C + T * cov * S, 1.0 - T * (1.0 - cov));
}
