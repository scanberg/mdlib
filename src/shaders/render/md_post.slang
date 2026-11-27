// md_post.slang -- the kernels of md_post.c: the display transform (viamd's
// tonemap/display.frag), the temporal resolve with the motion blur (its
// temporal.frag, velocity/blit_tilemax.frag and blit_neighbormax.frag) and
// the present with the sharpening and the marks (its sharpen and marks
// passes).
//
// Pixels are x right, y down, their centres at + 0.5 (md_render's). Layout
// rules of md_gpu.slang: no 3-vectors, 4-vectors at 16-byte offsets,
// pointers and handles at 8. Mirrored in md_post.c.

#include "../md_gpu.slang"

struct Args {
    float4 bg;          // xyz: the background, display-referred linear; w: the exposure
    float4 tm;          // x: the curve, y: gamma, z: feedback_min, w: feedback_max
    float4 pa, pb, pc;  // the projection: P00, P11, P20, P21; P30, P31, P22, P32; P23, P33, a pixel's size
                        // in view space at linear depth 1 (or at any: orthographic), 1 orthographic
    float4 jit;         // xy: this frame's jitter, md_render's (pixels; the image moved by it)
    float4 vis0, vis1;  // the marks: the selection's and the highlight's tint in view
    float4 hid0, hid1;  // fill where hidden
    float4 out0, out1;  // outline
    float4 g0, g1, g2, g3, g4, g5, g6;  // the highlight's groups 1..7 (rgb)
    float4 mk;          // x: how much of the highlight, y: the outline's width (output pixels),
                        // z: the dither's seed, w: 1 with marks
    uint   in_w, in_h;  // the frame's size (md_render's; the masks' in present)
    uint   out_w, out_h;// the output's
    uint   has_prev;    // resolve: a history to blend with
    uint   img_w, img_h;// present: the image's size
    uint   _p0;
    DescriptorHandle<Texture2D<float4>>   src;      // display: the scene; resolve: the frame; present: the image
    DescriptorHandle<Texture2D<float4>>   prev;     // resolve: the previous history
    DescriptorHandle<Texture2D<float>>    depth;    // resolve: md_render's, reverse Z, 0 where nothing
    DescriptorHandle<Texture2D<float2>>   motion;   // resolve: md_render's, pixels
    DescriptorHandle<Texture2D<float4>>   mask;     // present: md_render's masks
    DescriptorHandle<SamplerState>        bilinear; // clamp to edge
    DescriptorHandle<RWTexture2D<float4>> dst;
    float  sharpen;     // present: the unsharp mask's weight (0: none)
    uint   copy_shown;  // resolve, without the blur: the image shown (dst2) is the temporal result
    float4 mb;          // the motion blur: x its scale (0: none), y the noise's seed, z, w the tiles across, down
    float2* tmax;       // vel_tilemax's: a tile's fastest motion (frame pixels)
    float2* tiles;      // vel_neighbormax's: the fastest about a tile, which the resolve blurs along
    DescriptorHandle<RWTexture2D<float4>> dst2;     // resolve: the image shown, with the motion blur
};
MD_SHADER_ARGS(Args);

// ---- Display transform ----------------------------------------------------------

static const float WHITE_POINT       = 8.0;     // three stops above reference white
static const float HABLE_SCALE       = 2.0;     // Hable's published input scale and white point
static const float HABLE_WHITE_POINT = 11.2;

float3 neutral(float3 color) {
    const float start_compression = 0.8 - 0.04;
    const float desaturation = 0.15;
    float x = min(color.r, min(color.g, color.b));
    float offset = x < 0.08 ? x - 6.25 * x * x : 0.04;
    color -= offset;
    float peak = max(color.r, max(color.g, color.b));
    if (peak < start_compression) return color;
    const float d = 1.0 - start_compression;
    float new_peak = 1.0 - d * d / (peak + d - start_compression);
    color *= new_peak / peak;
    float g = 1.0 - 1.0 / (desaturation * (peak - new_peak) + 1.0);
    return lerp(color, float3(new_peak), g);
}

// Stephen Hill's fit (BakingLab's ACES.hlsl): sRGB => XYZ => D65_2_D60 =>
// AP1 => RRT_SAT, the RRT and ODT, then ODT_SAT => XYZ => D60_2_D65 => sRGB
float3 aces(float3 color) {
    float3 v = float3(dot(color, float3(0.59719, 0.35458, 0.04823)),
                      dot(color, float3(0.07600, 0.90834, 0.01566)),
                      dot(color, float3(0.02840, 0.13383, 0.83777)));
    float3 a = v * (v + 0.0245786) - 0.000090537;
    float3 b = v * (0.983729 * v + 0.4329510) + 0.238081;
    float3 r = a / b;
    return float3(dot(r, float3( 1.60475, -0.53108, -0.07367)),
                  dot(r, float3(-0.10208,  1.10813, -0.00605)),
                  dot(r, float3(-0.00327, -0.07276,  1.07602)));
}

float3 hable(float3 x) {
    const float A = 0.15, B = 0.50, C = 0.10, D = 0.20, E = 0.02, F = 0.30;
    return ((x * (A * x + C * B) + D * E) / (x * (A * x + B) + D * F)) - E / F;
}

float3 linear_to_srgb(float3 c) {
    c = saturate(c);
    return lerp(c * 12.92, 1.055 * pow(c, float3(1.0 / 2.4)) - 0.055, step(float3(0.0031308), c));
}

float3 curve(uint mode, float3 x) {
    if (mode == 1u) return neutral(x) / neutral(float3(WHITE_POINT)).r;
    if (mode == 2u) return aces(x) / aces(float3(WHITE_POINT)).r;
    if (mode == 3u) return hable(HABLE_SCALE * x) / hable(float3(HABLE_WHITE_POINT)).r;
    return x;
}

[shader("compute")]
[numthreads(8, 8, 1)]
void display(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    if (tid.x >= a.in_w || tid.y >= a.in_h) return;
    int2   q  = int2(tid.xy);
    float4 s  = a.src.Load(int3(q, 0));
    float  al = saturate(s.a);
    float3 scene = al > 1.0e-4 ? curve(uint(a.tm.x), max(s.rgb * (a.bg.w / al), float3(0.0))) : float3(0.0);
    float3 color = saturate(scene) * al + a.bg.xyz * (1.0 - al);
    a.dst[q] = float4(pow(linear_to_srgb(color), float3(2.2 / a.tm.y)), 1.0);
}

// ---- Temporal resolve -------------------------------------------------------------
// Playdead's TemporalReprojection (MIT) as viamd had it: the velocity of the
// nearest of the 3x3, the 3x3's colour box (rounded: the mean of the 3x3's and
// the cross's), the history clipped toward its centre, the feedback from the
// unbiased difference in luminance (Lottes); the colour samples unjittered.

static const float Z_BG = 1.0e9;
static const float FLT_EPS = 0.00000001;

float linear_depth(Args a, float d) {
    if (d == 0.0) return Z_BG;
    return -(a.pb.w - d * a.pc.y) / (d * a.pc.x - a.pb.z);
}

float luminance(float3 c) { return dot(c, float3(0.2125, 0.7154, 0.0721)); }

int2 clamp_in(Args a, int2 q) { return clamp(q, int2(0), int2(int(a.in_w), int(a.in_h)) - 1); }

float4 clip_aabb(float3 aabb_min, float3 aabb_max, float4 p, float4 q) {
    // clips toward the box's centre only (but fast)
    float3 p_clip = 0.5 * (aabb_max + aabb_min);
    float3 e_clip = 0.5 * (aabb_max - aabb_min) + FLT_EPS;
    float4 v_clip = q - float4(p_clip, p.w);
    float3 v_unit = v_clip.xyz / e_clip;
    float3 a_unit = abs(v_unit);
    float  ma = max(a_unit.x, max(a_unit.y, a_unit.z));
    return ma > 1.0 ? float4(p_clip, p.w) + v_clip / ma : q;
}

// PREFILTER. Geometry finer than a pixel comes and goes from one frame's
// samples to the next as the jitter moves them, and the history, clipped to
// what this frame caught, flickers. Where the 3x3 shows such detail the frame
// is taken from all nine samples, weighted by where each lies in the scene.
// How far c stands out from both a and b the same way: a peak or a dip one
// sample wide (0 on a step and on a slope)
float ridge(float a, float b, float c) {
    float da = a - c, db = b - c;
    return (da * db > 0.0) ? min(abs(da), abs(db)) : 0.0;
}

float ridge_both(int i, int j, int k, float l[9], float d[9], float2 unit) {
    return min(ridge(l[i], l[j], l[k]) * unit.x, ridge(d[i], d[j], d[k]) * unit.y);
}

// How much the 3x3 holds detail finer than a pixel, 0 to 1
float fine_detail(Args a, float l[9], float d[9]) {
    float  px   = a.pc.w > 0.5 ? a.pc.z : a.pc.z * d[4];
    float2 unit = float2(1.0 / 0.1, 1.0 / max(0.5 * px, 1e-8));
    float r = ridge_both(3, 5, 4, l, d, unit);
    r = max(r, ridge_both(1, 7, 4, l, d, unit));
    r = max(r, ridge_both(0, 8, 4, l, d, unit));
    r = max(r, ridge_both(2, 6, 4, l, d, unit));
    r = max(r, ridge_both(0, 2, 1, l, d, unit));
    r = max(r, ridge_both(6, 8, 7, l, d, unit));
    r = max(r, ridge_both(0, 6, 3, l, d, unit));
    r = max(r, ridge_both(2, 8, 5, l, d, unit));
    return smoothstep(0.5, 1.0, r);
}

// The nine samples weighted, a row and a column of weights, normalised
float4 weighted_3x3(float4 c[9], float3 wx, float3 wy) {
    float4 row0 = wx.x * c[0] + wx.y * c[1] + wx.z * c[2];
    float4 row1 = wx.x * c[3] + wx.y * c[4] + wx.z * c[5];
    float4 row2 = wx.x * c[6] + wx.y * c[7] + wx.z * c[8];
    return (wy.x * row0 + wy.y * row1 + wy.z * row2) / (dot(wx, float3(1.0)) * dot(wy, float3(1.0)));
}

static const float PREFILTER_K = 1.2;   // exp(-K d^2), d in the frame's pixels

// ---- Motion blur ------------------------------------------------------------------
// viamd's (temporal.frag's USE_MOTION_BLUR, after McGuire et al. 2012, "A
// Reconstruction Filter for Plausible Motion Blur"): the frame's motion
// reduced to tiles of VEL_TILE pixels, each the fastest of its nearest
// surface (vel_tilemax), then to the fastest about each tile, VEL_EXTENT tiles
// a side, of what is not behind it (vel_neighbormax): so a moving surface
// smears over what is behind it, not the other way. The resolve samples the
// frame along that motion, scaled, and blends it with the temporal result by
// how fast it is (2 to 15 output pixels) into the image shown; the history
// stays without it.

static const int   VEL_TILE   = 8;
static const int   VEL_EXTENT = 3;
static const int   BLUR_TAPS  = 5;      // a side
static const float TRUST_FULL = 2.0;    // output pixels of motion below which there is no blur
static const float TRUST_NONE = 15.0;   // and above which it is all blur

// A workgroup a tile (VEL_TILE = 8: a thread a pixel): the nearest depth
// over the tile, then over the pixels within 1 % of it the fastest motion,
// each a reduction in groupshared memory
groupshared float  vt_z[64];
groupshared float4 vt_v[64];        // xy: a motion, z: its length squared

[shader("compute")]
[numthreads(8, 8, 1)]
void vel_tilemax(uint3 tid : SV_DispatchThreadID, uint3 gid : SV_GroupID, uint gi : SV_GroupIndex) {
    Args a  = MD_ARGS;
    uint tw = uint(a.mb.z), th = uint(a.mb.w);
    bool inside = tid.x < a.in_w && tid.y < a.in_h;
    float z = inside ? linear_depth(a, a.depth.Load(int3(int2(tid.xy), 0))) : 1e30;
    vt_z[gi] = z;
    GroupMemoryBarrierWithGroupSync();
    [ForceUnroll] for (uint r = 32u; r > 0u; r >>= 1u) {
        if (gi < r) vt_z[gi] = min(vt_z[gi], vt_z[gi + r]);
        GroupMemoryBarrierWithGroupSync();
    }
    float zmin = vt_z[0];
    float eps  = 0.01 * zmin;       // a thickness, so the surface is not dropped for its quantisation
    float2 v   = inside && z <= zmin + eps ? a.motion.Load(int3(int2(tid.xy), 0)) : float2(0.0);
    vt_v[gi] = float4(v, dot(v, v), 0.0);
    GroupMemoryBarrierWithGroupSync();
    [ForceUnroll] for (uint r = 32u; r > 0u; r >>= 1u) {
        if (gi < r && vt_v[gi + r].z > vt_v[gi].z) vt_v[gi] = vt_v[gi + r];
        GroupMemoryBarrierWithGroupSync();
    }
    if (gi == 0u && gid.x < tw && gid.y < th) a.tmax[gid.y * tw + gid.x] = vt_v[0].xy;
}

// The linear depth at a tile's centre pixel
float tile_depth(Args a, int2 t) {
    return linear_depth(a, a.depth.Load(int3(clamp_in(a, t * VEL_TILE + VEL_TILE / 2), 0)));
}

[shader("compute")]
[numthreads(8, 8, 1)]
void vel_neighbormax(uint3 tid : SV_DispatchThreadID) {
    Args a  = MD_ARGS;
    int  tw = int(a.mb.z), th = int(a.mb.w);
    if (int(tid.x) >= tw || int(tid.y) >= th) return;
    int2   t  = int2(tid.xy);
    float  zc = tile_depth(a, t);
    float2 mv = float2(0.0);
    float  mv2 = 0.0;
    for (int j = -VEL_EXTENT; j <= VEL_EXTENT; ++j) {
        for (int i = -VEL_EXTENT; i <= VEL_EXTENT; ++i) {
            int2 n = clamp(t + int2(i, j), int2(0), int2(tw - 1, th - 1));
            if (tile_depth(a, n) > zc) continue;
            float2 v = a.tmax[n.y * tw + n.x];
            if (dot(v, v) > mv2) { mv = v; mv2 = dot(v, v); }
        }
    }
    a.tiles[t.y * tw + t.x] = mv;
}

// The tiles' motion at a point of the frame (pixels), bilinear between the
// tiles' centres
float2 tiles_at(Args a, float2 px) {
    int    tw = int(a.mb.z), th = int(a.mb.w);
    float2 tp = px / float(VEL_TILE) - 0.5;
    float2 b  = floor(tp);
    float2 f  = tp - b;
    int2   t0 = clamp(int2(b), int2(0), int2(tw - 1, th - 1));
    int2   t1 = clamp(int2(b) + 1, int2(0), int2(tw - 1, th - 1));
    float2 v00 = a.tiles[t0.y * tw + t0.x], v10 = a.tiles[t0.y * tw + t1.x];
    float2 v01 = a.tiles[t1.y * tw + t0.x], v11 = a.tiles[t1.y * tw + t1.x];
    return lerp(lerp(v00, v10, f.x), lerp(v01, v11, f.x), f.y);
}

// in [-1, 1] (Playdead's PDsrand)
float pd_srand(float2 n) { return frac(sin(dot(n, float2(12.9898, 78.233))) * 43758.5453) * 2.0 - 1.0; }

// The frame along vt (uv) about uv: 2 BLUR_TAPS + 1 samples over half of it
// each way, weighted by a triangle, offset by a noise that breaks up the
// trailing shells into grain
float4 sample_motion(Args a, float2 uv, float2 vt) {
    float2 v    = 0.5 * vt;
    float  sr   = pd_srand(uv + float2(a.mb.y, a.mb.y));
    float2 vtap = v / float(BLUR_TAPS);
    float2 pos0 = uv + sr * vtap * 0.5;
    float4 acc  = float4(0.0);
    float  wsum = 0.0;
    for (int i = -BLUR_TAPS; i <= BLUR_TAPS; ++i) {
        float w = float(BLUR_TAPS - abs(i) + 1);
        acc  += w * a.src.SampleLevel(a.bilinear, pos0 + float(i) * vtap, 0.0);
        wsum += w;
    }
    return acc / wsum;
}

// The image shown at the output pixel tid from the temporal result c: c, or
// blended with the frame along the motion
void write_shown(Args a, uint3 tid, float2 g_px, float2 j, float4 c) {
    if (a.mb.x <= 0.0) {
        if (a.copy_shown != 0u) a.dst2[int2(tid.xy)] = c;
        return;
    }
    float2 in_size  = float2(float(a.in_w), float(a.in_h));
    float2 out_size = float2(float(a.out_w), float(a.out_h));
    float2 sv    = a.mb.x * tiles_at(a, g_px);                 // frame pixels
    float  speed = length(sv * out_size / in_size);            // output pixels
    float  trust = 1.0 - clamp(speed - TRUST_FULL, 0.0, TRUST_NONE - TRUST_FULL) / (TRUST_NONE - TRUST_FULL);
    float4 cm    = trust < 1.0 ? sample_motion(a, (g_px - j) / in_size, sv / in_size) : c;
    a.dst2[int2(tid.xy)] = lerp(cm, c, trust);
}

[shader("compute")]
[numthreads(8, 8, 1)]
void resolve(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    if (tid.x >= a.out_w || tid.y >= a.out_h) return;
    float2 in_size  = float2(float(a.in_w), float(a.in_h));
    float2 out_size = float2(float(a.out_w), float(a.out_h));
    bool   upscale  = a.in_w != a.out_w || a.in_h != a.out_h;
    float2 tc   = (float2(tid.xy) + 0.5) / out_size;
    // texel q holds the scene at q + 0.5 + j: the image moved by the jitter
    float2 j    = -a.jit.xy;
    float2 g_px = tc * in_size;
    float2 g_px_cur = g_px, g_offset = float2(0.0);
    if (upscale) {
        // the texel whose sample lies nearest the output pixel's centre
        float2 q0 = floor(g_px - j);
        g_px_cur = q0 + 0.5;
        g_offset = q0 + 0.5 + j - g_px;
    }
    int2 qc = int2(floor(g_px_cur));

    // the 3x3: depths, colours, row by row from (-1, -1)
    float  d[9];
    float4 c[9];
    int    imin = 0;
    [ForceUnroll] for (int i = 0; i < 9; ++i) {
        int2 q = clamp_in(a, qc + int2(i % 3 - 1, i / 3 - 1));
        d[i] = linear_depth(a, a.depth.Load(int3(q, 0)));
        c[i] = a.src.Load(int3(q, 0));
        if (d[i] < d[imin]) imin = i;
    }
    // the velocity of the nearest, as uv of the output
    float2 m  = a.motion.Load(int3(clamp_in(a, qc + int2(imin % 3 - 1, imin / 3 - 1)), 0));
    float2 vel = m / in_size;

    // this frame at the output pixel's centre, unjittered
    float4 texel0 = a.src.SampleLevel(a.bilinear, (g_px - j) / in_size, 0.0);
    float  confidence = 1.0;
    float  l[9];
    [ForceUnroll] for (int i = 0; i < 9; ++i) l[i] = luminance(c[i].rgb);
    float  k = fine_detail(a, l, d);
    if (upscale) {
        // the samples weighted by their distance from the output pixel's
        // centre in output pixels (Blackman-Harris as a Gaussian); how near a
        // sample came is how much this frame counts
        float2 scale = out_size / in_size;
        float3 ox = g_offset.x + float3(-1.0, 0.0, 1.0);
        float3 oy = g_offset.y + float3(-1.0, 0.0, 1.0);
        float3 dx = ox * scale.x, dy = oy * scale.y;
        float3 wx = exp(-2.29 * dx * dx), wy = exp(-2.29 * dy * dy);
        texel0     = weighted_3x3(c, wx, wy);
        confidence = min(dot(wx, float3(1.0)) * dot(wy, float3(1.0)), 1.0);
        if (k > 0.0) {
            texel0     = lerp(texel0, weighted_3x3(c, exp(-PREFILTER_K * ox * ox), exp(-PREFILTER_K * oy * oy)), k);
            confidence = lerp(confidence, 1.0, k);
        }
    } else if (k > 0.0) {
        float3 dx = float3(-1.0, 0.0, 1.0) + j.x, dy = float3(-1.0, 0.0, 1.0) + j.y;
        texel0 = lerp(texel0, weighted_3x3(c, exp(-PREFILTER_K * dx * dx), exp(-PREFILTER_K * dy * dy)), k);
    }
    if (a.has_prev == 0u) {
        a.dst[int2(tid.xy)] = texel0;
        write_shown(a, tid, g_px, j, texel0);
        return;
    }
    float4 texel1 = a.prev.SampleLevel(a.bilinear, tc - vel, 0.0);

    // the colour box of the 3x3, rounded
    float4 cmin = c[0], cmax = c[0], cavg = c[0];
    [ForceUnroll] for (int i = 1; i < 9; ++i) { cmin = min(cmin, c[i]); cmax = max(cmax, c[i]); cavg += c[i]; }
    cavg /= 9.0;
    float4 cmin5 = min(c[1], min(c[3], min(c[4], min(c[5], c[7]))));
    float4 cmax5 = max(c[1], max(c[3], max(c[4], max(c[5], c[7]))));
    float4 cavg5 = (c[1] + c[3] + c[4] + c[5] + c[7]) / 5.0;
    cmin = 0.5 * (cmin + cmin5);
    cmax = 0.5 * (cmax + cmax5);
    cavg = 0.5 * (cavg + cavg5);
    texel1 = clip_aabb(cmin.xyz, cmax.xyz, clamp(cavg, cmin, cmax), texel1);

    float lum0 = luminance(texel0.rgb), lum1 = luminance(texel1.rgb);
    float unbiased = 1.0 - abs(lum0 - lum1) / max(lum0, max(lum1, 0.2));
    float fb = lerp(a.tm.z, a.tm.w, unbiased * unbiased);
    float4 res = upscale ? lerp(texel1, texel0, (1.0 - fb) * confidence) : lerp(texel0, texel1, fb);
    a.dst[int2(tid.xy)] = res;
    write_shown(a, tid, g_px, j, res);
}

// ---- Present ------------------------------------------------------------------------

float hash(float2 p) {
    float3 p3 = frac(float3(p.xyx) * 0.1031);
    p3 += dot(p3, p3.yzx + 33.33);
    return frac((p3.x + p3.y) * p3.z);
}

// A texel of the masks: what of each channel is set (1 or 0)
float4 mask_set(Args a, int2 q) {
    float4 m = a.mask.Load(int3(clamp_in(a, q), 0));
    return select(m > float4(0.0), float4(1.0), float4(0.0));
}

// The masks' set channels at a point (mask pixels), filtered bilinearly
float4 mask_at(Args a, float2 mp) {
    float2 b = floor(mp - 0.5);
    float2 f = mp - 0.5 - b;
    int2   q = int2(b);
    return lerp(lerp(mask_set(a, q), mask_set(a, q + int2(1, 0)), f.x),
                lerp(mask_set(a, q + int2(0, 1)), mask_set(a, q + int2(1, 1)), f.x), f.y);
}

// The mean of the four about the output pixel px, r output pixels out diagonally
float4 marks_about(Args a, float2 px, float r, float2 to_mask) {
    return 0.25 * (mask_at(a, (px + float2(-r, -r)) * to_mask) + mask_at(a, (px + float2(r, -r)) * to_mask) +
                   mask_at(a, (px + float2(-r,  r)) * to_mask) + mask_at(a, (px + float2(r,  r)) * to_mask));
}

// A mask texel's highlight tag: of the surface in view, else of anywhere; -1 none
int tag_of(float4 m) {
    if (m.z > 0.0) return int(round(m.z * 15.0)) - 8;
    if (m.w > 0.0) return int(round(m.w * 15.0)) - 8;
    return -1;
}

// The highlight's group at a point: its texel's, else of the 3x3 about it
// (the tent reaches that far) in view first
int group_about(Args a, float2 mp) {
    int2 q = clamp_in(a, int2(floor(mp)));
    int  t = tag_of(a.mask.Load(int3(q, 0)));
    if (t < 0) {
        int best = -1;
        for (int i = 0; i < 9; ++i) {
            float4 m = a.mask.Load(int3(clamp_in(a, q + int2(i % 3 - 1, i / 3 - 1)), 0));
            int    s = tag_of(m);
            if (s >= 0) best = max(best, s + (m.z > 0.0 ? 8 : 0));
        }
        t = best >= 0 ? best % 8 : -1;
    }
    return t;
}

float3 group_color(Args a, int g) {
    switch (g) {
    case 1: return a.g0.rgb;
    case 2: return a.g1.rgb;
    case 3: return a.g2.rgb;
    case 4: return a.g3.rgb;
    case 5: return a.g4.rgb;
    case 6: return a.g5.rgb;
    default: return a.g6.rgb;
    }
}

// Whether a mark is within reach of the group's 8 x 8 output pixels: in the
// mask texels its taps can read (marks_about's, at 0.5 and the outline's
// width, a texel more about them), 2 x 2 to a bilinear fetch at their shared
// corner (a set channel is 8 / 15 or more, so a quarter of it shows)
groupshared uint pr_marks;
bool marks_in_reach(Args a, uint2 gid, uint gi, float2 to_mask) {
    if (gi == 0u) pr_marks = 0u;
    GroupMemoryBarrierWithGroupSync();
    float  r   = max(0.5, a.mk.y - 0.5);
    float2 g0  = float2(gid * 8u);
    int2   top = int2(int(a.in_w), int(a.in_h)) - 1;
    int2   lo  = clamp(int2(floor((g0 + 0.5 - r) * to_mask - 0.5)) - 1, int2(0), top);
    int2   hi  = clamp(int2(floor((g0 + 7.5 + r) * to_mask - 0.5)) + 2, int2(0), top);
    int2   nb  = (hi - lo) / 2 + 1;                                     // the 2 x 2 blocks across, down
    float2 in_size = float2(float(a.in_w), float(a.in_h));
    bool   hit = false;
    for (int k = int(gi); k < nb.x * nb.y && !hit; k += 64) {
        float2 corner = float2(lo + 2 * int2(k % nb.x, k / nb.x) + 1);
        hit = any(a.mask.SampleLevel(a.bilinear, corner / in_size, 0.0) > float4(0.0));
    }
    if (hit) InterlockedOr(pr_marks, 1u);
    GroupMemoryBarrierWithGroupSync();
    return pr_marks != 0u;
}

[shader("compute")]
[numthreads(8, 8, 1)]
void present(uint3 tid : SV_DispatchThreadID, uint3 gid : SV_GroupID, uint gi : SV_GroupIndex) {
    Args a = MD_ARGS;
    float2 out_size = float2(float(a.out_w), float(a.out_h));
    float2 to_mask  = float2(float(a.in_w), float(a.in_h)) / out_size;
    bool   marks    = false;
    if (a.mk.w > 0.5) marks = marks_in_reach(a, gid.xy, gi, to_mask);     // the whole group, before any leaves
    if (tid.x >= a.out_w || tid.y >= a.out_h) return;
    float2 px = float2(tid.xy) + 0.5;
    bool   same = a.img_w == a.out_w && a.img_h == a.out_h;
    float4 c = same ? a.src.Load(int3(int2(tid.xy), 0)) : a.src.SampleLevel(a.bilinear, px / out_size, 0.0);
    // the sharpening (viamd's sharpen pass): an unsharp mask of the four
    // neighbours, (1 + w) c - w / 4 (l + r + t + b), at the output's size
    if (a.sharpen > 0.0 && same) {
        int2   q  = int2(tid.xy);
        int2   hi = int2(int(a.img_w), int(a.img_h)) - 1;
        float3 n  = a.src.Load(int3(clamp(q + int2(-1, 0), int2(0), hi), 0)).rgb + a.src.Load(int3(clamp(q + int2(1, 0), int2(0), hi), 0)).rgb +
                    a.src.Load(int3(clamp(q + int2(0, -1), int2(0), hi), 0)).rgb + a.src.Load(int3(clamp(q + int2(0, 1), int2(0), hi), 0)).rgb;
        c.rgb = (1.0 + a.sharpen) * c.rgb - 0.25 * a.sharpen * n;
    }
    if (marks) {
        float4 m    = marks_about(a, px, 0.5, to_mask);
        float4 ring = a.mk.y > 1.0 ? 0.5 * (m + marks_about(a, px, a.mk.y - 0.5, to_mask)) : m;
        if (max(max(ring.y, ring.w), max(m.y, m.w)) > 0.0) {
            float2 vis    = m.xz;
            float2 hidden = max(m.yw - m.xz, float2(0.0));
            float2 edge   = saturate(2.0 - 2.0 * abs(2.0 * ring.yw - 1.0));
            float3 v0 = a.vis0.rgb, h0 = a.hid0.rgb, o0 = a.out0.rgb;
            float3 v1 = a.vis1.rgb, h1 = a.hid1.rgb, o1 = a.out1.rgb;
            if (max(ring.w, m.w) > 0.0) {
                int g = group_about(a, px * to_mask);
                if (g > 0) { float3 gc = group_color(a, g); v1 = gc; h1 = gc; o1 = gc; }
            }
            c.rgb = lerp(c.rgb, v0, a.vis0.a * vis.x);
            c.rgb = lerp(c.rgb, h0, a.hid0.a * hidden.x);
            c.rgb = lerp(c.rgb, o0, a.out0.a * edge.x);
            float k = a.mk.x;
            c.rgb = lerp(c.rgb, v1, a.vis1.a * vis.y * k);
            c.rgb = lerp(c.rgb, h1, a.hid1.a * hidden.y * k);
            c.rgb = lerp(c.rgb, o1, a.out1.a * edge.y * k);
        }
    }
    float2 p = px + a.mk.z * 61.0;
    float  d = hash(p) + hash(p + float2(17.31, 41.17)) - 1.0;   // triangular in [-1, 1]
    a.dst[int2(tid.xy)] = float4(saturate(c.rgb + d / 255.0), 1.0);
}
