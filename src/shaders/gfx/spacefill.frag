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

void write_vis(ivec2 coord, float depth, uint payload) {
	ivec2 c = clamp(coord, ivec2(0,0), ivec2(ubo.render_width-1, ubo.render_height-1));
	uint idx = coord.y * ubo.render_width + coord.x;
	uint64_t cur_val = vis_buf[idx];
	float cur_depth = uintBitsToFloat(uint(cur_val >> 32));
	if (depth < cur_depth) {
		uint64_t val = (uint64_t(floatBitsToUint(depth)) << 32) | payload;
		atomicMin(vis_buf[idx], val);
	}
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

    float depth = (-ubo.view_to_clip[2][2] - ubo.view_to_clip[3][2] / (view_coord.z)) * 0.5 + 0.5;
    write_vis(ivec2(gl_FragCoord.xy), depth, IN.index);
}