#version 460

//#extension GL_NV_gpu_shader5 : enable
#extension GL_ARB_gpu_shader_int64 : enable
#extension GL_NV_shader_atomic_int64 : enable

#include "common.h"

#ifndef ORTHO
#define ORTHO 0
#endif

layout (binding = 0, std140) uniform UboBuffer {
    UniformData ubo;
};

layout(location = 0) in GS_FS {
    flat   vec4 view_sphere;
    smooth vec3 view_coord;
    flat   uint index;
#if ORTHO
    smooth vec2 uv;
#endif
} IN;

layout (binding = VIS_BINDING, std430) buffer VisBuffer {
	uint64_t vis_buf[];
};

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

    float depth = (-ubo.view_to_clip[2][2] - ubo.view_to_clip[3][2] / (view_coord.z)) * 0.5 + 0.5;

    uint depth_uint = floatBitsToUint(depth);
    uint payload = IN.index;
    uint64_t vis_data = (uint64_t(depth_uint) << 32) | payload;
    uint64_t vis_value = atomicMin(vis_buf[uint(gl_FragCoord.y) * ubo.render_width + uint(gl_FragCoord.x)], vis_data);
}