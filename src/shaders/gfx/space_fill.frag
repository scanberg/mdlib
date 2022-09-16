#version 460
#extension GL_EXT_conservative_depth : enable

struct UniformData {
    mat4 world_to_view;
    mat4 world_to_view_normal;
    mat4 world_to_clip;
    mat4 view_to_clip;
    mat4 view_to_world;
    mat4 clip_to_view;
    mat4 prev_world_to_clip;
    mat4 curr_view_to_prev_clip;
    vec4 jitter_uv;
    uint atom_mask;
    uint ortho;
    vec2 fle;
    float radius_scale;
};

struct VertexData {
    smooth vec3 view_coord;
    flat   vec4 view_sphere;
    flat   vec3 view_velocity;
    flat   vec4 color;
    flat   uint atom_idx;
    smooth vec2 uv;
};

layout(location = 0) out VertexData in_frag;

#ifdef GL_EXT_conservative_depth
layout (depth_greater) out float gl_FragDepth;
#endif

layout(location = 0) out vec4 out_color;
layout(location = 1) out vec3 out_normal;
layout(location = 2) out vec3 out_velocity;
layout(location = 3) out uint out_atom_index;

void main() {
    vec3 center = in_frag.view_sphere.xyz;
    float radius = in_frag.view_sphere.w;

if (u_ortho == 1U) {
    vec2 uv = in_frag.uv;
    float len = length(uv);
    if (len > 1.0) discard;
    vec3 view_coord = in_frag.view_coord.xyz + vec3(0, 0, radius * (sqrt(1.0 - len*len)));
} else {
    vec3 m = -center;
    vec3 d = normalize(in_frag.view_coord);
    float r = radius;
    float b = dot(m, d);
    float c = dot(m, m) - r*r;
    float discr = b*b - c;
    if (discr < 0.0) discard;
    float t = -b -sqrt(discr);
    vec3 view_coord = d * t;
}

    vec4 clip_coord = u_view_to_clip * vec4(view_coord, 1);
    gl_FragDepth = (clip_coord.z / clip_coord.w) * 0.5 + 0.5;

    vec3 view_normal = (view_coord - center) / radius;
    vec4 color = in_frag.color;
    vec3 view_velocity = in_frag.view_velocity;
    uint atom_index = in_frag.atom_idx;

    out_color = color;
    out_normal = view_normal;
    out_velocity = view_velocity;
    out_atom_index = atom_index;
}