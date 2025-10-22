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
    uint u_atom_base_index;
    uint u_bond_base_index;
    uint _pad0;
    float u_radius;
    float u_max_d2;
};

in Fragment {
    flat vec3  capsule_pa;
    flat vec3  capsule_pb;
    flat float capsule_ra;
    flat vec4  color[2];
    flat uint  atom_picking_idx[2];
    flat uint  bond_picking_idx;
    smooth vec3 view_vel;
    smooth vec3 view_pos;
} in_frag;

// Source from Ingo Quilez (https://www.shadertoy.com/view/Xt3SzX)
// Returns the ray scalar 't'
float capIntersect( in vec3 ro, in vec3 rd, in vec3 pa, in vec3 pb, in float r )
{
    vec3  ba = pb - pa;
    vec3  oa = ro - pa;

    float baba = dot(ba,ba);
    float bard = dot(ba,rd);
    float baoa = dot(ba,oa);
    float rdoa = dot(rd,oa);
    float oaoa = dot(oa,oa);

    float a = baba      - bard*bard;
    float b = baba*rdoa - baoa*bard;
    float c = baba*oaoa - baoa*baoa - r*r*baba;
    float h = b*b - a*c;
    if( h>=0.0 )
    {
        float t = (-b-sqrt(h))/a;
        float y = baoa + t*bard;
        // body
        if( y>0.0 && y<baba ) return t;
        // caps
        vec3 oc = (y<=0.0) ? oa : ro - pb;
        b = dot(rd,oc);
        c = dot(oc,oc) - r*r;
        h = b*b - c;
        if( h>0.0 ) return -b - sqrt(h);
    }
    return -1.0;
}

// compute normal and segment t value (how far along the)
vec4 capNormalAndSeg( in vec3 pos, in vec3 a, in vec3 b, in float r )
{
    vec3  ba = b - a;
    vec3  pa = pos - a;
    float h = clamp(dot(pa,ba)/dot(ba,ba),0.0,1.0);
    vec3  n = (pa - h*ba)/r;
    return vec4(n, h);
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
    vec3  pa = in_frag.capsule_pa;
    vec3  pb = in_frag.capsule_pb;
    float r  = in_frag.capsule_ra;

    float t = capIntersect(ro, rd, pa, pb, r);
    if (t < 0.0) {
        discard;
        return;
    }
    vec3 pos = ro + rd * t;
    vec4 normal_seg = capNormalAndSeg(pos, pa, pb, r);

    vec3 view_normal = normal_seg.xyz;
    float seg_t = normal_seg.w;

    int side = int(seg_t + 0.5);
    vec3 view_coord = rd * t;
    vec4 clip_coord = u_view_to_clip * vec4(view_coord, 1);
    gl_FragDepth = (clip_coord.z / clip_coord.w) * 0.5 + 0.5;

    vec4 color = in_frag.color[side];
    vec3 view_velocity = in_frag.view_vel;
    uint picking_index = abs(0.5 - seg_t) > 0.25 ? in_frag.atom_picking_idx[side] : in_frag.bond_picking_idx;

    write_fragment(view_coord, view_velocity, view_normal, color, picking_index);
}