// md_compose.slang -- the compose of md_compose.c: md_render's transparent
// layers shaded (viamd's compose_deferred.frag, both models) and composited
// front to back. Layout rules of md_gpu.slang: no 3-vectors, 4-vectors at
// 16-byte offsets, pointers at 8. Mirrored in md_compose.c.

#include "../md_gpu.slang"

struct Args {
    float4 ip0, ip1, ip2, ip3;  // the inverse projection's columns (the layout check models no arrays)
    float4 light;               // xyz: toward the key light (unit), w: its strength
    float4 dir_radiance;        // xyz, w: F0
    float4 sky;                 // xyz, w: roughness
    float4 ground;              // xyz, w: basic ambient
    float4 misc;                // x: basic specular, y: 1 basic, z: 1 linear colours, w: 1 orthographic
    uint*  layers;              // md_render_targets_t.layers: colour, normal, depth planes
    uint   count;               // layers a pixel (K)
    uint   stride;              // words of a layer (pixels)
    uint   width;
    uint   height;
    uint   _pad0, _pad1;
    DescriptorHandle<RWTexture2D<float4>> target;
};
MD_SHADER_ARGS(Args);

static const float PI          = 3.1415926535;
static const float ONE_OVER_PI = 1.0 / 3.1415926535;

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

// compose_deferred.frag's shade, without occlusion and shadow
float3 shade_pbr(Args a, float3 color, float3 V, float3 N) {
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
    Lo += (1.0 - Fe) * albedo * environment_irradiance(a, N) + Fe * environment_radiance(a, R, roughness);
    return Lo;
}

static const float BASIC_SHININESS = 48.0;

// compose_deferred.frag's shade_basic, without shadow
float3 shade_basic(Args a, float3 color, float3 V, float3 N) {
    float3 L = a.light.xyz;
    float  w = 0.5 + 0.5 * dot(N, L);
    float  light = a.light.w;
    float  amb = a.ground.w;
    float3 c = saturate(color) * (amb + light * (1.0 - amb) * w * w);
    c /= max(1.0, max(c.r, max(c.g, c.b)));
    float3 H = normalize(L + V);
    float  s = light * a.misc.x * pow(max(dot(N, H), 0.0), BASIC_SHININESS);
    c = lerp(c, float3(1.0), saturate(s));
    return a.misc.z > 0.5 ? c : srgb_to_linear(c);
}

float3 oct_decode(uint bits) {
    float2 e = float2(f16tof32(bits & 0xFFFFu), f16tof32(bits >> 16));
    float3 n = float3(e, 1.0 - abs(e.x) - abs(e.y));
    if (n.z < 0.0) {
        float2 s = select(n.xy >= 0.0, float2(1.0), float2(-1.0));
        n.xy = (1.0 - abs(n.yx)) * s;
    }
    return normalize(n);
}

[shader("compute")]
[numthreads(8, 8, 1)]
void compose_layers(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    if (tid.x >= a.width || tid.y >= a.height) return;
    uint   p = tid.y * a.width + tid.x;
    float3 C = float3(0.0);
    float  T = 1.0;
    // the pixel's direction toward the eye: through its centre, at the near plane and at depth 0.5
    float2 ndc = float2((float(tid.x) + 0.5) / float(a.width) * 2.0 - 1.0, 1.0 - (float(tid.y) + 0.5) / float(a.height) * 2.0);
    for (uint k = 0u; k < a.count; ++k) {
        float d = asfloat(a.layers[2u * a.count * a.stride + k * a.stride + p]);
        if (d == 0.0) break;
        float4 h = a.ip0 * ndc.x + a.ip1 * ndc.y + a.ip2 * d + a.ip3;
        float3 P = h.xyz / h.w;
        float3 V = a.misc.w > 0.5 ? float3(0.0, 0.0, 1.0) : -normalize(P);
        uint   c8 = a.layers[k * a.stride + p];
        float4 c  = float4(float(c8 & 0xFFu), float((c8 >> 8) & 0xFFu), float((c8 >> 16) & 0xFFu), float(c8 >> 24)) / 255.0;
        float3 N  = oct_decode(a.layers[a.count * a.stride + k * a.stride + p]);
        float3 L  = a.misc.y > 0.5 ? shade_basic(a, c.rgb, V, N) : shade_pbr(a, c.rgb, V, N);
        C += T * c.a * L;
        T *= 1.0 - c.a;
    }
    a.target[tid.xy] = float4(C, 1.0 - T);
}
