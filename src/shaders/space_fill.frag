#version 330 core
#extension GL_ARB_conservative_depth : enable

#ifndef ORTHO
#define ORTHO 0
#endif

in GS_FS {
    smooth vec3 view_coord;
    flat vec4 view_sphere;
    flat vec3 view_velocity;
    flat vec4 color;
    flat uint atom_idx;
    smooth vec2 uv;
} in_frag;

layout (std140) uniform ubo {
    mat4 u_world_to_view;
    mat4 u_world_to_view_normal;
    mat4 u_world_to_clip;
    mat4 u_view_to_clip;
    mat4 u_view_to_world;
    mat4 u_clip_to_view;
    mat4 u_prev_world_to_clip;
    mat4 u_curr_view_to_prev_clip;
    vec4 u_jitter_uv;
    uint u_atom_mask;
    uint _pad0;
    uint _pad1;
    uint _pad2;
    float u_radius_scale;
};

#ifdef GL_EXT_conservative_depth
layout (depth_greater) out float gl_FragDepth;
#endif

#pragma EXTRA_SRC

void main() {
    vec3 center = in_frag.view_sphere.xyz;
    float radius = in_frag.view_sphere.w;

#if ORTHO
    vec2 uv = in_frag.uv;
    float len = length(uv);
    if (len > 1.0) discard;
    vec3 view_coord = in_frag.view_coord.xyz + vec3(0, 0, radius * (sqrt(1.0 - len*len)));
#else
    vec3 m = -center;
    vec3 d = normalize(in_frag.view_coord);
    float r = radius;
    float b = dot(m, d);
    float c = dot(m, m) - r*r;
    float discr = b*b - c;
    if (discr < 0.0) discard;
    float t = -b -sqrt(discr);
    vec3 view_coord = d * t;
#endif

    vec4 clip_coord = u_view_to_clip * vec4(view_coord, 1);
    gl_FragDepth = (clip_coord.z / clip_coord.w) * 0.5 + 0.5;

    vec3 view_normal = (view_coord - center) / radius;
    vec4 color = in_frag.color;
    vec3 view_velocity = in_frag.view_velocity;
    uint atom_index = in_frag.atom_idx;

    write_fragment(view_coord, view_velocity, view_normal, color, atom_index);
}