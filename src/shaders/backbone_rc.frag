#version 330 core
#extension GL_ARB_conservative_depth : enable

// Segments of the backbone with the stadium profile sphere traced within the quads of backbone_rc.vert (see there for
// the solid). The distance to it is Lindow's: from the point p, the closest point pa on the core of ring 0 and pb on
// that of ring 1, then the closest point q on the segment between them, and |p - q| less the radius blended along it.
// That is exact for a capsule (no flat sides) and close for a band, which twists little over a segment.
//
// The normal is p - q, that of the capsule which holds the hit, with the component along the tangent of the spline
// (blended between the rings) taken out: the segments are straight, and their own normals would crease the surface where
// they meet (a helix bends by more than 15 degrees a segment). Taken perpendicular to the tangent of the curve, the normal
// is continuous over the joints, as that of the mesh is. Not at the round end of a chain, which keeps its own.
//
// The march starts where the ray enters the box of the segment and ends where it leaves it (see main for how). A ray
// that starts inside the segment is a miss.

#ifndef ORTHO
#define ORTHO 0
#endif

#define MAX_STEPS 64
#define MAX_REFINE 10

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

    vec4  u_scale;
    uvec4 u_res;
    uvec4 u_rings;
};

in Fragment {
    flat vec3  core[4];
    flat vec2  radius;
    flat vec3  vel[2];
    flat vec3  tangent[2];
    flat uvec2 open_end;
    flat vec3  box_lo;
    flat vec3  box_hi;
    flat vec4  color;
    flat uint  picking_idx;
    smooth vec3 view_pos;
} in_frag;

// The quad of backbone_rc.vert is in front of the segment, so the depth only increases
#ifdef GL_ARB_conservative_depth
layout (depth_greater) out float gl_FragDepth;
#endif

#pragma EXTRA_SRC

// The closest point to p on the segment [a, b], and its parameter
vec3 closest_on_segment(vec3 p, vec3 a, vec3 b, out float h) {
    vec3 ba = b - a;
    float bb = dot(ba, ba);
    h = bb > 0.0 ? clamp(dot(p - a, ba) / bb, 0.0, 1.0) : 0.0;
    return a + ba * h;
}

// The distance from p to the segment, the closest point on the axis of its capsule and the parameter along it
float segment_distance(vec3 p, out vec3 q, out float h) {
    float ha, hb;
    vec3 pa = closest_on_segment(p, in_frag.core[0], in_frag.core[1], ha);
    vec3 pb = closest_on_segment(p, in_frag.core[2], in_frag.core[3], hb);
    q = closest_on_segment(p, pa, pb, h);
    return length(p - q) - mix(in_frag.radius.x, in_frag.radius.y, h);
}

void main() {
    vec3 ro = in_frag.view_pos;
#if ORTHO
    vec3 rd = vec3(0, 0, -1);
#else
    vec3 rd = normalize(in_frag.view_pos);
#endif

    // The ray within the box
    vec3 inv_rd = 1.0 / rd;
    vec3 t_a = (in_frag.box_lo - ro) * inv_rd;
    vec3 t_b = (in_frag.box_hi - ro) * inv_rd;
    vec3 t_near = min(t_a, t_b);
    vec3 t_far  = max(t_a, t_b);
    float t0 = max(max(max(t_near.x, t_near.y), t_near.z), 0.0);
    float t1 = min(min(t_far.x, t_far.y), t_far.z);
    if (t0 > t1) {
        discard;
    }

    // Sphere trace until the distance changes sign, then find the root in that bracket (Illinois regula falsi). A plain
    // sphere tracer that stops where the distance is below a threshold converges slowly at grazing angles, where the ray
    // runs along a flat face: it ran out of steps there (holes) or stopped anywhere along the stretch within the
    // threshold (a depth off by threshold / sin(angle), which made the overlapping segments fight). With the sign change
    // the hit is the surface itself, and a step is never less than h_min, which bounds the steps through the box.
    // The distance is 1-Lipschitz (its gradient is that of the distance to the point q), so a step by it never crosses
    // the surface; the steps by h_min may cross pieces thinner than that along the ray, which are the slivers at the
    // silhouettes of the round edges (a percent of the radius).
    float r_min = min(in_frag.radius.x, in_frag.radius.y);
    float h_min = max(0.25 * r_min, length(in_frag.box_hi - in_frag.box_lo) / float(MAX_STEPS - 16));
    float tol   = 1.0e-4 * r_min;

    vec3  q;
    float h;
    float t  = t0;
    float d  = segment_distance(ro + rd * t, q, h);
    if (d < 0.0) {
        // The ray starts inside (the eye, or the near plane cutting through the segment): seen from within it is open,
        // as a mesh is to back face culling
        discard;
    }

    bool  bracketed = false;
    float ta = t, fa = d, tb = t, fb = d;
    for (int i = 0; i < MAX_STEPS; ++i) {
        float tn = min(t + max(d, h_min), t1);
        float dn = segment_distance(ro + rd * tn, q, h);
        if (dn < 0.0) {
            ta = t;  fa = d;
            tb = tn; fb = dn;
            bracketed = true;
            break;
        }
        if (tn >= t1) break;
        t = tn;
        d = dn;
    }
    if (!bracketed) {
        discard;
    }

    // Illinois: regula falsi which halves the value kept at an end that stays twice, so it converges from both sides
    int side = 0;
    t = tb;
    for (int j = 0; j < MAX_REFINE; ++j) {
        t = (ta * fb - tb * fa) / (fb - fa);
        float ft = segment_distance(ro + rd * t, q, h);
        if (abs(ft) < tol) break;
        if (ft < 0.0) {
            tb = t; fb = ft;
            if (side == -1) fa *= 0.5;
            side = -1;
        } else {
            ta = t; fa = ft;
            if (side == 1) fb *= 0.5;
            side = 1;
        }
    }

    vec3 view_coord  = ro + rd * t;
    vec3 view_normal = normalize(view_coord - q);
    bool at_end = (h <= 0.0 && in_frag.open_end.x != 0u) || (h >= 1.0 && in_frag.open_end.y != 0u);
    if (!at_end) {
        vec3 tng = normalize(mix(in_frag.tangent[0], in_frag.tangent[1], h));
        vec3 n   = view_normal - tng * dot(view_normal, tng);
        float n2 = dot(n, n);
        if (n2 > 1.0e-6) view_normal = n * inversesqrt(n2);
    }
    vec3 view_vel    = mix(in_frag.vel[0], in_frag.vel[1], h);

    vec4 clip_coord = u_view_to_clip * vec4(view_coord, 1);
    gl_FragDepth = (clip_coord.z / clip_coord.w) * 0.5 + 0.5;

    write_fragment(view_coord, view_vel, view_normal, in_frag.color, in_frag.picking_idx);
}
