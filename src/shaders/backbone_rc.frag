#version 330 core
#extension GL_ARB_conservative_depth : enable

// Segments of the backbone with the stadium profile ray cast within the quads of backbone_rc.vert, in closed form.
//
// The solid Lindow et al. describe is the union of the capsules between corresponding points of the cores of the two
// rings: the ruled surface between the cores (a bilinear patch) thickened by the radius. The patch is a quadric, but
// its offset by r is not, and has no practical closed form. Its sphere tracing needs a distance to it, and Lindow's
// (to the segment between the closest points of the two cores) is not one: it can overestimate, so the march needs a
// least step, and misses slivers thinner than that at grazing angles.
//
// Here the segment is the union of pieces which each meet a ray in closed form:
//   - the flat faces: two bilinear patches, the core ends moved by +-r along the normal of the flat sides at their ring
//     (y = z x x, as the mesh of backbone.vert places its profile). Reshetov's ray / patch, a quadratic.
//   - the rounded sides: the round cones (convex hulls) of the spheres at the ends of the cores, (a0, r0) to (b0, r1)
//     and (a1, r0) to (b1, r1). Quilez's ray / round cone, a quadratic.
//   - where the band ends (open_end): at the end of a chain the rungs, the capsules along the cores, and the round cones
//     whole; where the neighbor is not drawn a flat cap, the plane of the ring within the stadium, the segment cut there
//     as elsewhere (a pass which draws the highlighted atoms only gives the same depth where it draws as the pass of
//     all, which drew the neighbor). Elsewhere the faces of the next segment meet these at the shared ring exactly: the
//     segments of a chain close each other.
// The side walls between the faces (the core ends moved by up to r) lie within the round cones, so the union is closed
// and the first hit is the nearest hit with any piece: no iteration. The segment is cut at the planes of its rings
// (through the middle of the core, along the tangent), but where the chain ends: a hit with a piece past one is
// dropped, the neighbor there holds the same part of the band (the round cones of the two meet in nearly the same
// ellipse, a few thousandths of r apart). So the segments of a chain do not overlap, and the quads of backbone_rc.vert
// bound their own parts only. The faces lie between the planes, the round cones reach past them. Against the thickened ruled patch (dense sampling)
// it is off by 0.01 Angstrom at full detail and a few hundredths at the coarsest (the march: up to 0.13).
// A round tube (l = 0 at both rings) is the round cone alone, exactly.
//
// The normal is that of the piece, turned by the rotation which takes the axis of the segment to the tangent of the
// spline there (blended between the rings): the segments are straight, and their own normals would crease the surface
// where they meet. Turned with the curve, the normal is continuous over the joints, as that of the mesh is, and keeps
// the tilt the twist of a band or a change of its profile gives it (taking out its component along the tangent would
// lose that). Not past the round end of a band, which keeps its own.

#ifndef ORTHO
#define ORTHO 0
#endif

#define NO_HIT 3.4e38

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
    flat vec3  normal[2];
    flat vec3  vel[2];
    flat vec3  tangent[2];
    flat uvec2 open_end;        // 0 not, 1 round (a chain), 2 flat (a neighbor not drawn)
    flat vec4  color;
    flat uint  picking_idx;
    smooth vec3 view_pos;
} in_frag;

// The quad of backbone_rc.vert is in front of the segment, so the depth only increases
#ifdef GL_ARB_conservative_depth
layout (depth_greater) out float gl_FragDepth;
#endif

#pragma EXTRA_SRC

// Ray (ro, rd) against the bilinear patch q00, q10, q01, q11 (u from q00 to q10, v from q00 to q01): the nearest hit
// with u and v in [0, 1] and 0 < t < t_hit, and the normal of the patch there (either side). Reshetov, "Cool Patches: A
// Geometric Approach to Ray/Bilinear Patch Intersections" (Ray Tracing Gems, 2019): the ray meets the line of the patch
// at u where the two are coplanar, a quadratic in u, then t and v on that line.
void patch_hit(vec3 ro, vec3 rd, vec3 q00, vec3 q10, vec3 q01, vec3 q11, inout float t_hit, inout vec3 n_hit, inout vec3 n_side) {
    vec3 e10 = q10 - q00;
    vec3 e11 = q11 - q10;
    vec3 e00 = q01 - q00;
    vec3 p00 = q00 - ro;
    vec3 p10 = q10 - ro;
    float a = dot(cross(p00, rd), e00);
    float c = dot(cross(e10, q01 - q11), rd);
    float b = dot(cross(p10, rd), e11) - a - c;
    float det = b * b - 4.0 * a * c;
    if (det < 0.0) return;
    det = sqrt(det);
    float u1, u2;
    if (c == 0.0) {
        // A trapezoid: one root
        u1 = -a / b;
        u2 = -1.0;
    } else {
        // The roots without cancellation: q / c and a / q
        float q = -0.5 * (b + (b < 0.0 ? -det : det));
        u1 = q / c;
        u2 = a / q;
    }
    for (int i = 0; i < 2; ++i) {
        float u = i == 0 ? u1 : u2;
        if (!(u >= 0.0 && u <= 1.0)) continue;
        vec3 pa = mix(p00, p10, u);     // The line of the patch at u: pa + v pb
        vec3 pb = mix(e00, e11, u);
        vec3 n  = cross(rd, pb);
        float dd = dot(n, n);
        if (!(dd > 0.0)) continue;
        n = cross(n, pa);
        float t = dot(n, pb) / dd;
        float v = dot(n, rd) / dd;
        if (t > 0.0 && t < t_hit && v >= 0.0 && v <= 1.0) {
            t_hit  = t;
            n_hit  = cross(mix(e10, q11 - q01, v), pb);
            n_side = n_hit;
        }
    }
}

// Whether a point is within the segment: between the planes of its rings, but where the chain ends
bool kept(vec3 p) {
    float tol = 1.0e-3 * (in_frag.radius.x + in_frag.radius.y);
    vec3  m0  = 0.5 * (in_frag.core[0] + in_frag.core[1]);
    vec3  m1  = 0.5 * (in_frag.core[2] + in_frag.core[3]);
    return (in_frag.open_end.x == 1u || dot(p - m0, in_frag.tangent[0]) >= -tol) &&
           (in_frag.open_end.y == 1u || dot(p - m1, in_frag.tangent[1]) <=  tol);
}

// Ray against the flat cap at a ring: its plane (through the middle of the core [a, b], normal z), within r of the core
void cap_hit(vec3 ro, vec3 rd, vec3 a, vec3 b, float r, vec3 z, inout float t_hit, inout vec3 n_hit, inout vec3 n_side) {
    float dz = dot(rd, z);
    if (abs(dz) < 1.0e-8) return;
    vec3  m = 0.5 * (a + b);
    float t = dot(m - ro, z) / dz;
    if (!(t > 0.0 && t < t_hit)) return;
    vec3  p  = ro + rd * t;
    vec3  ba = b - a;
    float h  = clamp(dot(p - a, ba) / max(dot(ba, ba), 1.0e-12), 0.0, 1.0);
    vec3  q  = p - (a + ba * h);
    // Within r, not the ellipses the round cones of the sides make cut at a slant (wider by a few thousandths of r): a
    // cap must never reach out of the surface, which would change the depth where the pass of all drew the segment
    if (dot(q, q) <= r * r && kept(p)) {
        t_hit  = t;
        n_hit  = z;
        n_side = z;
    }
}

// Ray (ro, rd) against the round cone of the spheres (pa, ra) and (pb, rb), their convex hull (a capsule for ra = rb):
// the nearest hit with 0 < t < t_hit, within the segment (kept), and the outward normal there. Quilez,
// https://iquilezles.org/articles/intersectors. Its nearest hit past a plane of the segment drops it: the ray entered
// the band there through the neighbor.
void round_cone_hit(vec3 ro, vec3 rd, vec3 pa, vec3 pb, float ra, float rb, inout float t_hit, inout vec3 n_hit, inout vec3 n_side) {
    vec3  ba = pb - pa;
    vec3  oa = ro - pa;
    vec3  ob = ro - pb;
    float rr = ra - rb;
    float m0 = dot(ba, ba);
    float m1 = dot(ba, oa);
    float m2 = dot(ba, rd);
    float m3 = dot(rd, oa);
    float m5 = dot(oa, oa);
    float m6 = dot(ob, rd);
    float m7 = dot(ob, ob);
    float d2 = m0 - rr * rr;
    // Otherwise one sphere holds the other, or they coincide: the spheres alone
    if (d2 > 0.0) {
        float k2 = d2 - m2 * m2;
        float k1 = d2 * m3 - m1 * m2 + m2 * rr * ra;
        float k0 = d2 * m5 - m1 * m1 + m1 * rr * ra * 2.0 - m0 * ra * ra;
        float h  = k1 * k1 - k0 * k2;
        // The spheres are within the cone: a ray which misses it misses them
        if (h < 0.0) return;
        float t = (-sqrt(h) - k1) / k2;
        float y = m1 - ra * rr + t * m2;
        if (y > 0.0 && y < d2) {
            if (t > 0.0 && t < t_hit && kept(ro + rd * t)) {
                t_hit  = t;
                n_hit  = d2 * (oa + t * rd) - ba * y;
                n_side = n_hit;
            }
            return;
        }
    }
    // The spheres: the nearer hit of the two
    float h1 = m3 * m3 - m5 + ra * ra;
    float h2 = m6 * m6 - m7 + rb * rb;
    float t1 = h1 > 0.0 ? -m3 - sqrt(h1) : NO_HIT;
    float t2 = h2 > 0.0 ? -m6 - sqrt(h2) : NO_HIT;
    t1 = t1 > 0.0 ? t1 : NO_HIT;
    t2 = t2 > 0.0 ? t2 : NO_HIT;
    float t = min(t1, t2);
    if (t < t_hit && kept(ro + rd * t)) {
        t_hit = t;
        n_hit = t1 <= t2 ? oa + t * rd : ob + t * rd;
        // That of the cone through it: where a segment meets the next, its sphere is cut at the plane of the ring, and
        // turned with the curve the normal of the sphere would turn twice
        float y = m1 - ra * rr + t * m2;
        n_side = d2 > 0.0 ? d2 * (oa + t * rd) - ba * y : n_hit;
    }
}

void main() {
    vec3 ro = in_frag.view_pos;
#if ORTHO
    vec3 rd = vec3(0, 0, -1);
#else
    vec3 rd = normalize(in_frag.view_pos);
#endif

    vec3  a0 = in_frag.core[0];
    vec3  a1 = in_frag.core[1];
    vec3  b0 = in_frag.core[2];
    vec3  b1 = in_frag.core[3];
    float r0 = in_frag.radius.x;
    float r1 = in_frag.radius.y;

    float t  = NO_HIT;
    vec3  n  = vec3(0.0);   // Of the surface
    vec3  ns = vec3(0.0);   // Of the band along the segment, to be turned with the curve

    // A side, or the whole of a round tube
    round_cone_hit(ro, rd, a0, b0, r0, r1, t, n, ns);

    vec3 ca = a1 - a0;
    vec3 cb = b1 - b0;
    float rs = r0 + r1;
    if (dot(ca, ca) + dot(cb, cb) > 1.0e-8 * rs * rs) {
        round_cone_hit(ro, rd, a1, b1, r0, r1, t, n, ns);

        // The faces towards the ray: a ray enters the solid through a face only against its normal. Both where the
        // normals of the rings disagree about it (grazing).
        vec3  y0 = in_frag.normal[0] * r0;
        vec3  y1 = in_frag.normal[1] * r1;
        float s0 = dot(rd, y0);
        float s1 = dot(rd, y1);
        if (s0 <= 0.0 || s1 <= 0.0) patch_hit(ro, rd, a0 + y0, a1 + y0, b0 + y1, b1 + y1, t, n, ns);
        if (s0 >= 0.0 || s1 >= 0.0) patch_hit(ro, rd, a0 - y0, a1 - y0, b0 - y1, b1 - y1, t, n, ns);

        if (in_frag.open_end.x == 1u) round_cone_hit(ro, rd, a0, a1, r0, r0, t, n, ns);
        if (in_frag.open_end.y == 1u) round_cone_hit(ro, rd, b0, b1, r1, r1, t, n, ns);
    }

    // The flat caps, last: a cap which is the nearest hit keeps its own normal
    float t_before = t;
    if (in_frag.open_end.x == 2u) cap_hit(ro, rd, a0, a1, r0, -in_frag.tangent[0], t, n, ns);
    if (in_frag.open_end.y == 2u) cap_hit(ro, rd, b0, b1, r1,  in_frag.tangent[1], t, n, ns);
    bool cap = t < t_before;

    if (t >= NO_HIT) {
        discard;
    }

    vec3 view_coord  = ro + rd * t;

    // Along the segment, from the middle of the core of ring 0 to that of ring 1
    vec3  m0 = 0.5 * (a0 + a1);
    vec3  ax = 0.5 * (b0 + b1) - m0;
    float h  = dot(view_coord - m0, ax) / max(dot(ax, ax), 1.0e-12);
    bool at_end = cap || (h < 0.0 && in_frag.open_end.x == 1u) || (h > 1.0 && in_frag.open_end.y == 1u);
    h = clamp(h, 0.0, 1.0);
    vec3 view_normal = normalize(at_end ? n : ns);
    if (dot(view_normal, rd) > 0.0) view_normal = -view_normal;
    if (!at_end && dot(ax, ax) > 1.0e-12) {
        // Rodrigues, without the angle: v + w x v + w x (w x v) / (1 + a.b), w = a x b
        vec3  a = normalize(ax);
        vec3  b = normalize(mix(in_frag.tangent[0], in_frag.tangent[1], h));
        vec3  w = cross(a, b);
        float cab = dot(a, b);
        if (cab > -0.9) {
            vec3 wv = cross(w, view_normal);
            view_normal = normalize(view_normal + wv + cross(w, wv) / (1.0 + cab));
        }
    }
    vec3 view_vel = mix(in_frag.vel[0], in_frag.vel[1], h);

    vec4 clip_coord = u_view_to_clip * vec4(view_coord, 1);
    gl_FragDepth = (clip_coord.z / clip_coord.w) * 0.5 + 0.5;

    write_fragment(view_coord, view_vel, view_normal, in_frag.color, in_frag.picking_idx);
}
