#version 330 core
#extension GL_ARB_conservative_depth : enable

#ifndef ORTHO
#define ORTHO 0
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
    float u_radius;
};

in Fragment {
    flat uint picking_idx[2];
    flat vec4 color[2];
    flat vec4 capsule_center_radius;
    flat vec4 capsule_axis_length;
    smooth vec3 view_vel;
    smooth vec3 view_pos;
} in_frag;

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

#ifdef GL_EXT_conservative_depth
layout (depth_greater) out float gl_FragDepth;
#endif

#pragma EXTRA_SRC

void main() {
#if ORTHO
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
    vec4 clip_coord = u_view_to_clip * vec4(view_coord, 1);
    gl_FragDepth = (clip_coord.z / clip_coord.w) * 0.5 + 0.5;

    vec4 color = in_frag.color[side];
    vec3 view_velocity = in_frag.view_vel;
    uint atom_index = in_frag.picking_idx[side];

    write_fragment(view_coord, view_velocity, view_normal, color, atom_index);
}