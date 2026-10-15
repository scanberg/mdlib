#version 330 core
#extension GL_ARB_conservative_depth : enable

// Slabs raycast within the quads of slab.vert: the Minkowski sum of a convex polygon (z = 0 in the frame of the slab)
// with a sphere of radius u_radius, its faces moved in to the half thickness u_plate. Everything is exact and closed
// form: the rim is a capsule along each edge, the faces are the planes z = +-u_plate within the polygon, and the nearest
// hit of the ray with any of them is the surface. A slab of n corners costs n capsule and two plane intersections.
//
// The polygon is taken to be convex (a ring of a sugar or a base, projected on its plane): a point is within it when it
// is to the left of every edge, which go counter clockwise. A concave one would get its plate filled in at the concave
// corners, under the rim.

#ifndef ORTHO
#define ORTHO 0
#endif

#define MAX_CORNERS 7

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
    float u_plate;
    float u_max_d2;
    uint  u_slab_base_index;
};

in Fragment {
    flat vec3  center;
    flat vec3  normal;
    flat vec3  tangent;
    flat vec4  corner[4];
    flat uint  count;
    flat vec4  color;
    flat uint  picking;
    flat vec3  view_vel;
    smooth vec3 view_pos;
} in_frag;

// The quad of slab.vert is in front of the slab, so the depth only increases
#ifdef GL_ARB_conservative_depth
layout (depth_greater) out float gl_FragDepth;
#endif

#pragma EXTRA_SRC

vec2 corner_pos(int k) {
    vec4 c = in_frag.corner[k / 2];
    return (k % 2) == 0 ? c.xy : c.zw;
}

float cross2(vec2 a, vec2 b) {
    return a.x * b.y - a.y * b.x;
}

// Ray (ro, rd) against the capsule of the segment [pa, pb] and radius r. rd is unit length. Returns t of the hit where
// the ray enters it, or -1 for none. From Inigo Quilez (https://iquilezles.org/articles/intersectors).
float capIntersect(vec3 ro, vec3 rd, vec3 pa, vec3 pb, float r) {
    vec3  ba   = pb - pa;
    vec3  oa   = ro - pa;
    float baba = dot(ba, ba);
    float bard = dot(ba, rd);
    float baoa = dot(ba, oa);
    float rdoa = dot(rd, oa);
    float oaoa = dot(oa, oa);
    float a = baba - bard * bard;
    float b = baba * rdoa - baoa * bard;
    float c = baba * oaoa - baoa * baoa - r * r * baba;
    float h = b * b - a * c;
    if (h >= 0.0) {
        float t = (-b - sqrt(h)) / a;
        float y = baoa + t * bard;
        // The body
        if (y > 0.0 && y < baba) return t;
        // The caps
        vec3 oc = (y <= 0.0) ? oa : ro - pb;
        b = dot(rd, oc);
        c = dot(oc, oc) - r * r;
        h = b * b - c;
        if (h > 0.0) return -b - sqrt(h);
    }
    return -1.0;
}

void main() {
    // The ray of the fragment, from its point on the quad, which is in front of the slab
    vec3 ro = in_frag.view_pos;
#if ORTHO
    vec3 rd = vec3(0, 0, -1);
#else
    vec3 rd = normalize(in_frag.view_pos);
#endif

    // Frame of the slab (local -> view)
    vec3 nrm   = in_frag.normal;
    vec3 tng   = in_frag.tangent;
    mat3 basis = mat3(tng, cross(nrm, tng), nrm);
    vec3 lo    = (ro - in_frag.center) * basis;     // transpose(basis) * v: view -> local
    vec3 ld    = rd * basis;

    int   n = int(in_frag.count);
    float r = u_radius;
    float h = clamp(u_plate, 0.0, r);

    float t_hit = 3.4e38;
    vec3  n_hit = vec3(0.0);

    // The plates: the face the ray comes from, within the polygon
    if (abs(ld.z) > 1.0e-6) {
        float side = lo.z > 0.0 ? 1.0 : -1.0;
        float t = (side * h - lo.z) / ld.z;
        if (t > 0.0) {
            vec2 pt = lo.xy + ld.xy * t;
            bool inside = true;
            for (int k = 0; k < MAX_CORNERS; ++k) {
                if (k >= n) break;
                vec2 a = corner_pos(k);
                vec2 b = corner_pos((k + 1) % n);
                if (cross2(b - a, pt - a) < 0.0) {
                    inside = false;
                }
            }
            if (inside) {
                t_hit = t;
                n_hit = vec3(0.0, 0.0, side);
            }
        }
    }

    // The rim
    for (int k = 0; k < MAX_CORNERS; ++k) {
        if (k >= n) break;
        vec3 a = vec3(corner_pos(k), 0.0);
        vec3 b = vec3(corner_pos((k + 1) % n), 0.0);
        float t = capIntersect(lo, ld, a, b, r);
        if (t > 0.0 && t < t_hit) {
            vec3  p  = lo + ld * t;
            vec3  ba = b - a;
            float s  = clamp(dot(p - a, ba) / max(dot(ba, ba), 1.0e-12), 0.0, 1.0);
            t_hit = t;
            n_hit = (p - (a + ba * s)) / r;
        }
    }

    if (t_hit >= 3.0e38) {
        discard;
    }

    vec3 lp          = lo + ld * t_hit;
    vec3 view_coord  = in_frag.center + basis * lp;
    vec3 view_normal = normalize(basis * n_hit);

    vec4 clip_coord = u_view_to_clip * vec4(view_coord, 1);
    gl_FragDepth = (clip_coord.z / clip_coord.w) * 0.5 + 0.5;

    write_fragment(view_coord, in_frag.view_vel, view_normal, in_frag.color, in_frag.picking);
}
