#version 460

#include "common.h"

#ifndef ORTHO
#define ORTHO 0
#endif

layout (binding = 0, std140) uniform UboBuffer {
    UniformData ubo;
};

layout(location = 0) in GS_FS {
    flat   vec4 view_sphere;
    flat   vec3 view_velocity;
    flat   uint color;
    flat   uint index;
    smooth vec3 view_coord;
#if ORTHO
    smooth vec2 uv;
#endif
} IN;

layout (depth_greater) out float gl_FragDepth;

layout(location = COLOR_TEX_ATTACHMENT)  out vec4 out_color;
layout(location = NORMAL_TEX_ATTACHMENT) out vec2 out_normal;
layout(location = SS_VEL_TEX_ATTACHMENT) out vec2 out_velocity;
layout(location = INDEX_TEX_ATTACHMENT)  out vec4 out_index;

vec2 encode_normal (vec3 n) {
   float p = sqrt(n.z * 8 + 8);
   return n.xy / p + 0.5;
}

vec4 encode_index(uint index) {
    return vec4(
        (index & 0x000000FFU) >> 0U,
        (index & 0x0000FF00U) >> 8U,
        (index & 0x00FF0000U) >> 16U,
        (index & 0xFF000000U) >> 24U) / 255.0;
}

void main() {
    vec3 center = IN.view_sphere.xyz;
    float radius = IN.view_sphere.w;
    vec3 view_coord;

#if ORTHO
    vec2 uv = IN.uv;
    float len = length(uv);
    if (len > 1.0) discard;
    view_coord = IN.view_coord.xyz + vec3(0, 0, radius * (sqrt(1.0 - len*len)));
#else
    vec3 m = -center;
    vec3 d = normalize(IN.view_coord);
    float r = radius;
    float b = dot(m, d);
    float c = dot(m, m) - r*r;
    float discr = b*b - c;
    if (discr < 0.0) discard;
    float t = -b -sqrt(discr);
    view_coord = d * t;
#endif

    vec4 clip_coord = ubo.view_to_clip * vec4(view_coord, 1);
    gl_FragDepth = (clip_coord.z / clip_coord.w) * 0.5 + 0.5;

    vec3 view_normal = (view_coord - center) / radius;
    vec4 color = unpackUnorm4x8(IN.color);
    vec3 view_velocity = IN.view_velocity;
    uint index = IN.index;

    out_color = color;
    out_normal = encode_normal(view_normal);
    out_velocity = vec2(0,0);
    out_index = encode_index(index);
}