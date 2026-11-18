// render_window.slang: md_render_window's one pass, a fullscreen triangle that
// shades a view's targets into the window's image. The look of
// md_render_snap: the base colour lit by its normal from over the camera's
// shoulder, the selection's visible part tinted and its outline (hidden too)
// drawn. A look at the targets, not viamd's shading.

#include "../../src/shaders/md_gpu.slang"

struct Args {
    DescriptorHandle<Texture2D<float4>> color;      // RGBA8, sRGB-encoded base colour; a: covered
    DescriptorHandle<Texture2D<float2>> normal;     // view-space normal, octahedral
    DescriptorHandle<Texture2D<float4>> mask;       // r: state 0 visible, g: state 0 anywhere
    uint width, height;
    uint linear_out;                                // the image is *_SRGB: write linear
    uint _pad;
};
MD_SHADER_ARGS(Args);

struct FullV {
    float4 pos : SV_Position;
};

[shader("vertex")]
FullV full_vs(uint vid : SV_VertexID) {
    FullV o;
    float2 p = float2(float((vid << 1) & 2u), float(vid & 2u));
    o.pos = float4(p * 2.0 - 1.0, 0.0, 1.0);
    if (MD_ARGS.width == 0u) o.pos = float4(0.0, 0.0, -2.0, 1.0);  // reads the arguments: both stages must declare the same struct
    return o;
}

float srgb_to_lin(float c) { return c <= 0.04045 ? c / 12.92 : pow((c + 0.055) / 1.055, 2.4); }
float lin_to_srgb(float c) { c = saturate(c); return c <= 0.0031308 ? c * 12.92 : 1.055 * pow(c, 1.0 / 2.4) - 0.055; }

float3 oct_decode(float2 e) {
    float3 n = float3(e, 1.0 - abs(e.x) - abs(e.y));
    if (n.z < 0.0) n.xy = (1.0 - abs(n.yx)) * float2(n.x >= 0.0 ? 1.0 : -1.0, n.y >= 0.0 ? 1.0 : -1.0);
    return normalize(n);
}

[shader("fragment")]
float4 compose_fs(FullV i) : SV_Target {
    Args  a  = MD_ARGS;
    uint2 px = uint2(i.pos.xy);
    float3 o = float3(0.13, 0.13, 0.15);                        // background, sRGB-encoded
    float4 c = a.color[px];
    if (c.a > 0.0) {
        const float3 L = normalize(float3(-0.35, 0.45, 0.82));  // view space, towards the viewer and up-left
        float3 n = oct_decode(a.normal[px]);
        float shade = 0.22 + 0.78 * max(dot(n, L), 0.0);
        o = float3(lin_to_srgb(srgb_to_lin(c.r) * shade), lin_to_srgb(srgb_to_lin(c.g) * shade), lin_to_srgb(srgb_to_lin(c.b) * shade));
    }
    float4 m = a.mask[px];
    if (m.r > 0.0) o = o * 0.6 + float3(0.95, 0.95, 0.1) * 0.4;  // the selection, visible part tinted
    if (m.g > 0.0) {                                            // its outline, hidden parts too
        bool edge = false;
        if (px.x + 1u < a.width  && a.mask[px + uint2(1, 0)].g == 0.0) edge = true;
        if (px.x > 0u            && a.mask[px - uint2(1, 0)].g == 0.0) edge = true;
        if (px.y + 1u < a.height && a.mask[px + uint2(0, 1)].g == 0.0) edge = true;
        if (px.y > 0u            && a.mask[px - uint2(0, 1)].g == 0.0) edge = true;
        if (edge) o = float3(1.0, 0.85, 0.1);
    }
    if (a.linear_out != 0u) o = float3(srgb_to_lin(o.r), srgb_to_lin(o.g), srgb_to_lin(o.b));
    return float4(o, 1.0);
}
