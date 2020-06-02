#version 330 core
#extension GL_ARB_conservative_depth : enable
#extension GL_ARB_shading_language_packing : enable

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
    uint  u_atom_mask;

    float u_radius;
};

in Fragment {
    flat vec4 view_velocity[2];
    flat uint picking_idx[2];
    flat vec4 color[2];
    flat vec4 capsule_center_radius;
    flat vec4 capsule_axis_length;
    smooth vec3 view_pos;
} in_frag;

#ifdef GL_EXT_conservative_depth
layout (depth_greater) out float gl_FragDepth;
#endif
layout(location = 0) out vec4 out_color;
layout(location = 1) out vec4 out_normal;
layout(location = 2) out vec4 out_ss_vel;
layout(location = 3) out vec4 out_picking;

// Source from Ingo Quilez (https://www.shadertoy.com/view/Xt3SzX)
float intersect_capsule(in vec3 ro, in vec3 rd, in vec3 cc, in vec3 ca, float cr,
                      float ch, out vec3 normal, out float seg_t)  // cc center, ca orientation axis, cr radius, ch height
{
    vec3 oc = ro - cc;
    ch *= 0.5;

    float card = dot(ca, rd);
    float caoc = dot(ca, oc);

    float a = 1.0 - card * card;
    float b = dot(oc, rd) - caoc * card;
    float c = dot(oc, oc) - caoc * caoc - cr * cr;
    float h = b * b - a * c;
    if (h < 0.0) return -1.0;
    float t = (-b - sqrt(h)) / a;

    float y = caoc + t * card;
    seg_t = clamp(y * 0.5 + 0.5, 0.0, 1.0);

    // body
    if (abs(y) < ch) {
        normal = normalize(oc + t * rd - ca * y);
        return t;
    }

    // caps
    float sy = sign(y);
    oc = ro - (cc + sy * ca * ch);
    b = dot(rd, oc);
    c = dot(oc, oc) - cr * cr;
    h = b * b - c;
    if (h > 0.0) {
        t = -b - sqrt(h);
        normal = normalize(oc + rd * t);
        return t;
    }

    return -1.0;
}

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
#ifdef ORTHO
    vec3 ro = vec3(in_frag.view_pos.xy, 0);
    vec3 rd = vec3(0,0,-1);
#else
    vec3 ro = vec3(0);
    vec3 rd = normalize(in_frag.view_pos);
#endif
    vec3  cc = in_frag.capsule_center_radius.xyz;
    float cr = in_frag.capsule_center_radius.w;
    vec3  ca = in_frag.capsule_axis_length.xyz;
    float ch = in_frag.capsule_axis_length.w;

    vec3 view_normal;
    float seg_t;
    float t = intersect_capsule(ro, rd, cc, ca, cr, ch, view_normal, seg_t);
    if (t < 0.0) {
        discard;
        return;
    }

    int side = int(seg_t + 0.5);

    vec3 view_coord = rd * t;
    vec4 view_velocity = mix(in_frag.view_velocity[0], in_frag.view_velocity[1], seg_t);
    vec4 color = in_frag.color[side];
    uint picking_idx = in_frag.picking_idx[side];

    vec4 curr_clip_coord = u_view_to_clip * vec4(view_coord, 1);

    vec3 prev_view_coord = view_coord - view_velocity.xyz;
    vec4 prev_clip_coord = u_curr_view_to_prev_clip * vec4(prev_view_coord, 1);

    // Remove jitter from samples to provide the actual velocity
    // This is crucial for the temporal reprojection to work properly
    // Otherwise the velocity will push the samples outside of the "reprojection" region
    vec2 curr_ndc = curr_clip_coord.xy / curr_clip_coord.w;
    vec2 prev_ndc = prev_clip_coord.xy / prev_clip_coord.w;
    vec2 ss_vel = (curr_ndc - prev_ndc) * 0.5 + (u_jitter_uv.xy - u_jitter_uv.zw);

    gl_FragDepth = (curr_clip_coord.z / curr_clip_coord.w) * 0.5 + 0.5;
    out_color   = color;
    out_normal  = encode_normal(view_normal);
    out_ss_vel  = vec4(ss_vel, 0, 0);
    out_picking = unpackUnorm4x8(picking_idx);
}