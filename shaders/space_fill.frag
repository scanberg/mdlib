#version 330 core
#extension GL_ARB_conservative_depth : enable
#extension GL_ARB_shading_language_packing : enable

#ifndef ATOM_VEL
#define ATOM_VEL 1
#endif

#ifndef ATOM_COL
#define ATOM_COL 1
#endif

#ifndef ATOM_IDX
#define ATOM_IDX 1
#endif

#ifndef VIEW_NORM
#define VIEW_NORM 1
#endif

in GS_FS {
    smooth vec3 view_coord;
    flat vec4 view_sphere;
#if ATOM_VEL
    flat vec3 view_velocity;
#endif
#if ATOM_COL
    flat vec4 color;
#endif
#if ATOM_IDX
    flat uint atom_idx;
#endif
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

    float u_radius_scale;
};

#ifdef GL_EXT_conservative_depth
layout (depth_greater) out float gl_FragDepth;
#endif
layout(location = 0) out vec4 out_color;
layout(location = 1) out vec4 out_normal;
layout(location = 2) out vec4 out_uv_vel;
layout(location = 3) out vec4 out_picking;

// https://aras-p.info/texts/CompactNormalStorage.html
vec4 encode_normal (vec3 n) {
    float p = sqrt(n.z*8+8);
    return vec4(n.xy/p + 0.5,0,0);
}

#ifndef GL_ARB_shading_language_packing
vec4 unpackUnorm4x8(in uint data) {
    return vec4(
        (data & 0x000000FFU) >> 0U,
        (data & 0x0000FF00U) >> 8U,
        (data & 0x00FF0000U) >> 16U,
        (data & 0xFF000000U) >> 24U) / 255.0;
}
#endif

void main() {
    vec3 center = in_frag.view_sphere.xyz;
    float radius = in_frag.view_sphere.w;
    vec3 view_dir = normalize(in_frag.view_coord);

    vec3 m = -center;
    vec3 d = view_dir;
    float r = radius;
    float b = dot(m, d);
    float c = dot(m, m) - r*r;
    float discr = b*b - c;
    if (discr < 0.0) discard;
    float t = -b -sqrt(discr);

    vec3 view_coord = d * t;
    vec4 clip_coord = u_view_to_clip * vec4(view_coord, 1);

    gl_FragDepth = (clip_coord.z / clip_coord.w) * 0.5 + 0.5;

#if ATOM_VEL
    // @NOTE: This is a bit obscure, we need to get the fragment clip coordinate in the previous frame
    // Since we raycast implicit geometry, we cannot directly transform any vertex and resort to interpolation.
    // We instead use the surface hit point (view_coord) and subtract the view_velocity and transform this to the clip space of the previous frame
    // This will give us the fragment coordinate in the previous frame
    vec3 prev_view_coord = view_coord - in_frag.view_velocity;
    vec4 prev_clip_coord = u_curr_view_to_prev_clip * vec4(prev_view_coord, 1);

    // @NOTE: Remove jitter from samples to provide the actual velocity
    // This is crucial for the temporal reprojection to work properly
    // Otherwise the velocity will push the samples outside of the "reprojection" region
    vec2 curr_ndc = clip_coord.xy / clip_coord.w;
    vec2 prev_ndc = prev_clip_coord.xy / prev_clip_coord.w;
    vec2 uv_vel = (curr_ndc - prev_ndc) * 0.5 + (u_jitter_uv.xy - u_jitter_uv.zw);
    out_uv_vel = vec4(uv_vel, 0, 0);
#endif
#if VIEW_NORM
    vec3 view_normal = (view_coord - center) / radius;
    out_normal = encode_normal(view_normal);
#endif
#if ATOM_IDX
    out_picking = unpackUnorm4x8(in_frag.atom_idx);
#endif
#if ATOM_COL
    out_color = in_frag.color;
#endif
}