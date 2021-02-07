#version 330 core
#extension GL_ARB_shading_language_packing : enable

#ifndef ATOM_COL
#define ATOM_COL 1
#endif

#ifndef ATOM_VEL
#define ATOM_VEL 0
#endif

#ifndef ATOM_IDX
#define ATOM_IDX 1
#endif

#ifndef VIEW_NORM
#define VIEW_NORM 1
#endif

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
    float u_width_scale;
    float u_height_scale;
};


#if ATOM_COL || ATOM_VEL || ATOM_IDX || VIEW_NORM
in Fragment {
#if ATOM_VEL
    smooth vec4 curr_clip_coord;
    smooth vec4 prev_clip_coord;
#endif
#if ATOM_COL
    smooth vec4 color;
#endif
#if VIEW_NORM
    smooth vec3 view_normal;
#endif
#if ATOM_IDX
    flat   uint picking_idx;
#endif
} in_frag;
#endif

layout(location = 0) out vec4 out_color;
layout(location = 1) out vec4 out_normal;
layout(location = 2) out vec4 out_ss_vel;
layout(location = 3) out vec4 out_picking;

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
#if ATOM_COL
    out_color   = in_frag.color;
#endif
#if VIEW_NORM
    out_normal  = encode_normal(normalize(in_frag.view_normal.xyz));
#endif
#if ATOM_VEL
    vec2 curr_ndc = in_frag.curr_clip_coord.xy / in_frag.curr_clip_coord.w;
    vec2 prev_ndc = in_frag.prev_clip_coord.xy / in_frag.prev_clip_coord.w;
    vec2 ss_vel = (curr_ndc - prev_ndc) * 0.5 + (u_jitter_uv.xy - u_jitter_uv.zw);
    out_ss_vel  = vec4(ss_vel, 0, 0);
#endif
#if ATOM_IDX
    out_picking = unpackUnorm4x8(in_frag.picking_idx);
#endif
}