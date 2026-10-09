#version 330 core
#extension GL_ARB_conservative_depth : enable
#extension GL_ARB_shading_language_packing : enable

// Bonds raycast within the quads of bond.vert.
//   WEAK_BONDS 0: licorice, a capsule (the cylinder of the bond with a half sphere at each end).
//   WEAK_BONDS 1: weak bonds, a dashed cylinder: u_dash_count dashes with flat ends, the first starting at pa and the
//                 last ending at pb, each covering the fraction u_dash_fill of the period of the pattern.

#ifndef ORTHO
#define ORTHO 0
#endif

#ifndef WEAK_BONDS
#define WEAK_BONDS 0
#endif

#define MODE_NEAREST 0
#define MODE_SMOOTH  1
#define MODE_UNIFORM 2

#ifndef GL_ARB_shading_language_packing
vec4 unpackUnorm4x8(uint p) {
    uvec4 iv = uvec4(p & 0xFFU, (p >> 8U) & 0xFFU, (p >> 16U) & 0xFFU, (p >> 24U) & 0xFFU);
    return vec4(iv) * (1.0f / 255.0f);
}
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
    uint u_backbone_base_index;

    float u_radius;
    float u_max_d2;
    int   u_mode;
    float u_sharpness;
    uint  u_uniform_color;
    uint  u_dash_count;
    float u_dash_fill;
};

in Fragment {
    flat vec3  pa;
    flat vec3  pb;
    flat float radius;
    flat vec4  color[2];
    flat vec3  view_vel[2];
    flat uint  atom_picking_idx[2];
    flat uint  bond_picking_idx;
    smooth vec3 view_pos;
} in_frag;

// The quad of bond.vert is in front of the bond, so the depth only increases
#ifdef GL_ARB_conservative_depth
layout (depth_greater) out float gl_FragDepth;
#endif

#pragma EXTRA_SRC

#if WEAK_BONDS
// Ray (ro, rd) against the dashed cylinder from pa to pb of radius r. rd is unit length.
// Returns (t, normal), t < 0 for a miss, and the coordinate along the axis in [0, 1] in seg_t.
//
// The ray is within the infinite cylinder for t in [t_in, t_out] (intersected with the slab of the pattern), where
// its coordinate along the axis is linear in t. If it enters within a dash it hits the side of the cylinder (or the
// flat end of the bond), otherwise it enters in a gap and hits the end of the next dash in its direction along the
// axis, if it reaches it before it leaves the cylinder. Constant time for any number of dashes.
vec4 dashed_cylinder_intersect(vec3 ro, vec3 rd, vec3 pa, vec3 pb, float r, uint count, float fill, out float seg_t) {
    vec3  ba  = pb - pa;
    float len = length(ba);
    vec3  ax  = ba / len;
    vec3  oc  = ro - pa;

    float k  = dot(rd, ax);         // Rate of the axial coordinate along the ray
    float y0 = dot(oc, ax);         // Axial coordinate at t = 0
    vec3  q  = rd - ax * k;         // Components perpendicular to the axis
    vec3  w  = oc - ax * y0;

    float k2 = dot(q, q);
    float k1 = dot(q, w);
    float k0 = dot(w, w) - r * r;

    float t_in  = -3.4e38;
    float t_out =  3.4e38;
    bool  side  = true;             // t_in is on the side of the cylinder (otherwise on an end plane)

    if (k2 > 0.0) {
        float h = k1 * k1 - k2 * k0;
        if (h < 0.0) return vec4(-1.0);
        h = sqrt(h);
        t_in  = (-k1 - h) / k2;
        t_out = (-k1 + h) / k2;
    } else if (k0 > 0.0) {
        // Along the axis, outside of the radius
        return vec4(-1.0);
    }

    if (k != 0.0) {
        float t0 = (0.0 - y0) / k;
        float t1 = (len - y0) / k;
        float ts_in  = min(t0, t1);
        float ts_out = max(t0, t1);
        if (ts_in > t_in) {
            t_in = ts_in;
            side = false;
        }
        t_out = min(t_out, ts_out);
    } else if (y0 < 0.0 || len < y0) {
        return vec4(-1.0);
    }

    if (t_in > t_out) return vec4(-1.0);

    // Pattern: dash i covers [i * period, i * period + dash], the last one ends at len
    float n      = float(max(count, 1u));
    float f      = clamp(fill, 0.0, 1.0);
    float period = len / (n - 1.0 + f);
    float dash   = period * f;

    float y = clamp(y0 + k * t_in, 0.0, len);
    float i = min(floor(y / period), n - 1.0);
    float t = t_in;

    if (y - i * period > dash) {
        // In the gap after dash i: the next dash begins at (i + 1) * period, dash i ends at i * period + dash
        if (k == 0.0) return vec4(-1.0);
        y = k > 0.0 ? (i + 1.0) * period : i * period + dash;
        t = (y - y0) / k;
        if (t > t_out) return vec4(-1.0);
        side = false;
    }

    seg_t = y / len;

    vec3 normal;
    if (side) {
        vec3 p = oc + rd * t;
        normal = (p - ax * dot(p, ax)) / r;
    } else {
        normal = k > 0.0 ? -ax : ax;
    }
    return vec4(t, normal);
}
#else
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
#endif

void main() {
    // The ray of the fragment, from its point on the quad, which is in front of the bond
    vec3 ro = in_frag.view_pos;
#if ORTHO
    vec3 rd = vec3(0, 0, -1);
#else
    vec3 rd = normalize(in_frag.view_pos);
#endif
    vec3  pa = in_frag.pa;
    vec3  pb = in_frag.pb;
    float r  = in_frag.radius;

#if WEAK_BONDS
    float seg_t;
    vec4 hit = dashed_cylinder_intersect(ro, rd, pa, pb, r, u_dash_count, u_dash_fill, seg_t);
    float t = hit.x;
    if (t < 0.0) {
        discard;
    }
    vec3 view_normal = hit.yzw;
#else
    float t = capIntersect(ro, rd, pa, pb, r);
    if (t < 0.0) {
        discard;
    }
    vec4 normal_seg = capNormalAndSeg(ro + rd * t, pa, pb, r);
    vec3 view_normal = normal_seg.xyz;
    float seg_t = normal_seg.w;
#endif

    vec3 view_coord = ro + rd * t;
    vec4 clip_coord = u_view_to_clip * vec4(view_coord, 1);
    gl_FragDepth = (clip_coord.z / clip_coord.w) * 0.5 + 0.5;

    int side = int(seg_t + 0.5);
    vec4 color;
    if (u_mode == MODE_NEAREST) {
        color = in_frag.color[side];
    } else if (u_mode == MODE_SMOOTH) {
        float s = u_sharpness * u_sharpness;
        float k = mix(1.0, 32.0, s);
        float a = pow(seg_t, k);
        float b = pow(1.0 - seg_t, k);
        color = mix(in_frag.color[0], in_frag.color[1], a / (a + b));
    } else {
        color = unpackUnorm4x8(u_uniform_color);
    }

    vec3 view_velocity = mix(in_frag.view_vel[0], in_frag.view_vel[1], seg_t);

#if WEAK_BONDS
    // A weak bond is not a bond of the molecule: it picks as the atom it is closest to
    uint picking_index = in_frag.atom_picking_idx[side];
#else
    uint picking_index = abs(0.5 - seg_t) > 0.25 ? in_frag.atom_picking_idx[side] : in_frag.bond_picking_idx;
#endif

    write_fragment(view_coord, view_velocity, view_normal, color, picking_index);
}
