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
//     (y = z x x, as the mesh of backbone.vert places its profile). Reshetov's ray / patch, a quadratic. Each where the
//     ray can enter it, against its normal somewhere: across a twisted band the normal of a face turns away from those
//     of the rings by more than ten degrees (it is normal to the sides, which slant where the band twists), so the
//     normals of the rings do not tell (a ray at a grazing angle entered a face they said it left, and went through the
//     band). The normal of the patch is bilinear over it, within the cone of those at its corners: within radius.z
//     of the mean normal of the rings (backbone_rc.vert).
//   - the rounded sides: as the mesh makes them, the ruled surfaces between the circles of radius r0 and r1 about the
//     core ends in the planes of the rings, through the points at the same angle about them (the profile of the mesh at
//     its rings, its quads between them). A ray meets the line of the surface at an angle where a trigonometric
//     polynomial of degree 2 vanishes (a quartic): Newton's method, from where the ray enters the oblique cone of the
//     circles about the side in the planes normal to the mean tangent w (a quadratic), which differs from the surface by
//     a hundredth of r or two. Two steps reach its first root to a few ten thousandths of an Angstrom (three to a
//     millionth, against the roots bracketed densely). The faces meet the sides along their lines at the angles of Y, tangent to them, and the
//     sides of the next segment meet these in the same circles: the band is closed and smooth at the joints but for the
//     turn of the curve. A quadric for a side cannot do that: about the side (from a0 to b0), which slants from the
//     tangent where the band twists (by 6 and 14 degrees at the two ends of a segment of a helix), its extents across
//     the planes of the two rings keep the ratio of the cosines of those slants (2.5 %), a step at every joint and a
//     saw tooth along the edge, where it bulges over the face and sinks below it in turn.
//   - where the band ends (open_end): at the end of a chain the rungs, the capsules along the cores (a sphere for a round
//     tube), and the sides cut at the plane normal to w there; where the neighbor is not drawn a flat cap, the plane of
//     the ring within the stadium, the segment cut there as elsewhere (a pass which draws the highlighted atoms only
//     gives the same depth where it draws as the pass of all, which drew the neighbor). Elsewhere the faces of the next
//     segment meet these at the shared ring exactly: the segments of a chain close each other.
// The side walls between the faces (the core ends moved by up to r) lie within the sides, so the union is closed and the
// first hit is the nearest hit with any piece: no iteration but that within a side. The segment is cut at the planes
// of its rings (through the middle of the core, along the tangent), but where the chain ends: a hit with a piece past
// one is dropped, the neighbor there holds the same part of the band. So the segments of a chain do not overlap, and
// the quads of backbone_rc.vert bound their own parts only. The faces and the sides lie between the planes.
// A round tube (l = 0 at both rings) is the oblique cone alone: its side is its axis, which follows the tangent, so
// the cones of two segments meet the plane of their ring in the same circle to second order in the turn (0.2 % of r at
// 7.7 degrees), and the ruled surface would pinch the tube where its frame twists (as the mesh does, an artifact of it).
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
    flat vec3  radius;          // r0, r1, the margin of the faces (backbone_rc.vert)
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
    // Within r: the sides meet the plane of the ring in the circles of the stadium
    if (dot(q, q) <= r * r && kept(p)) {
        t_hit  = t;
        n_hit  = z;
        n_side = z;
    }
}

// Ray (ro, rd) against the round cone of the spheres (pa, ra) and (pb, rb), their convex hull (a capsule for ra = rb, a
// sphere for pa = pb): t where the ray enters it (NO_HIT for a miss) and the outward normal there. Quilez,
// https://iquilezles.org/articles/intersectors
float round_cone_enter(vec3 ro, vec3 rd, vec3 pa, vec3 pb, float ra, float rb, out vec3 n) {
    n = vec3(0.0);
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
        if (h < 0.0) return NO_HIT;
        float t = (-sqrt(h) - k1) / k2;
        float y = m1 - ra * rr + t * m2;
        if (y > 0.0 && y < d2) {
            n = d2 * (oa + t * rd) - ba * y;
            return t;
        }
    }
    float h1 = m3 * m3 - m5 + ra * ra;
    float h2 = m6 * m6 - m7 + rb * rb;
    float t1 = h1 > 0.0 ? -m3 - sqrt(h1) : NO_HIT;
    float t2 = h2 > 0.0 ? -m6 - sqrt(h2) : NO_HIT;
    if (t1 >= NO_HIT && t2 >= NO_HIT) return NO_HIT;
    float t = min(t1, t2);
    n = t1 <= t2 ? oa + t * rd : ob + t * rd;
    return t;
}

// The normal of the cone of the round cone at p, also off its body (where the ray enters it through a cut)
vec3 round_cone_side_normal(vec3 p, vec3 pa, vec3 pb, float ra, float rb) {
    vec3  ba = pb - pa;
    float rr = ra - rb;
    float d2 = dot(ba, ba) - rr * rr;
    if (!(d2 > 0.0)) return p - (ra >= rb ? pa : pb);
    float y = dot(ba, p - pa) - ra * rr;
    return d2 * (p - pa) - ba * y;
}

// Narrows the interval s of the ray to the side of the plane through m with the inward normal z
void clip_to_plane(vec3 ro, vec3 rd, vec3 m, vec3 z, float tol, inout vec2 s) {
    float d = dot(rd, z);
    float f = dot(ro - m, z) + tol;     // >= 0 inside
    if (abs(d) < 1.0e-12) {
        if (f < 0.0) s = vec2(NO_HIT, -NO_HIT);
        return;
    }
    float t = -f / d;
    if (d > 0.0) s.x = max(s.x, t);
    else         s.y = min(s.y, t);
}

// The interval of the ray between the planes of the rings, but where the chain ends (kept)
vec2 segment_interval(vec3 ro, vec3 rd) {
    float tol = 1.0e-3 * (in_frag.radius.x + in_frag.radius.y);
    vec2  s   = vec2(-NO_HIT, NO_HIT);
    if (in_frag.open_end.x != 1u) clip_to_plane(ro, rd, 0.5 * (in_frag.core[0] + in_frag.core[1]),  in_frag.tangent[0], tol, s);
    if (in_frag.open_end.y != 1u) clip_to_plane(ro, rd, 0.5 * (in_frag.core[2] + in_frag.core[3]), -in_frag.tangent[1], tol, s);
    return s;
}

// Ray against the round cone cut to the segment (a rung): where it enters the part of it between the planes of the rings,
// the intersection of the two intervals, if nearer than t_hit. Within the cone a ray enters the segment through a cut,
// and there the cut is the surface, with the normal of the cone (n_side).
void round_cone_hit(vec3 ro, vec3 rd, vec3 pa, vec3 pb, float ra, float rb, inout float t_hit, inout vec3 n_hit, inout vec3 n_side) {
    vec3  n_in;
    float t_in = round_cone_enter(ro, rd, pa, pb, ra, rb, n_in);
    if (t_in >= NO_HIT) return;
    // Where it leaves: where the ray back from beyond the cone enters it
    float far = t_in + length(pb - pa) + 2.0 * (ra + rb) + 1.0e-3;
    vec3  n_out;
    float t_back = round_cone_enter(ro + rd * far, -rd, pa, pb, ra, rb, n_out);
    float t_out  = t_back < NO_HIT ? far - t_back : t_in;

    vec2  s  = segment_interval(ro, rd);
    float t0 = max(t_in, s.x);
    float t1 = min(t_out, s.y);
    if (!(t0 <= t1 && t0 > 0.0 && t0 < t_hit)) return;
    vec3 p = ro + rd * t0;
    t_hit  = t0;
    n_side = round_cone_side_normal(p, pa, pb, ra, rb);
    n_hit  = t0 == t_in ? n_in : n_side;
}

// The direction the circles of a side are normal to: the mean tangent of the segment w, or the side itself where it
// slants from that by more than 84 degrees (the cone would be cut nearly along its axis). As backbone_rc.vert.
vec3 side_axis(vec3 d, vec3 w) {
    float ld = length(d);
    float lw = length(w);
    return dot(d, w) > 0.1 * ld * lw ? w / lw : d / max(ld, 1.0e-12);
}

// Ray (ro, rd) against the oblique cone of the circles about pa + s d (d = pb - pa) of radius ra + s (rb - ra), in the
// planes normal to w (d.w > 0), s unbounded: the interval of t within it (empty: x > y). The point of the axis in the
// plane of p is at s = (p - pa).w / d.w, and |p - pa - s d|^2 = (ra + s (rb - ra))^2 is a quadratic in t.
vec2 oblique_cone_interval(vec3 ro, vec3 rd, vec3 pa, vec3 pb, float ra, float rb, vec3 w) {
    const vec2 none = vec2(NO_HIT, -NO_HIT);
    vec3  d  = pb - pa;
    float dw = dot(d, w);
    float dr = rb - ra;
    vec3  oa = ro - pa;
    float so = dot(oa, w) / dw;
    float sd = dot(rd, w) / dw;
    vec3  q0 = oa - d * so;         // p - pa - s d = q0 + t q1
    vec3  q1 = rd - d * sd;
    float R0 = ra + dr * so;        // The radius there, R0 + t R1
    float R1 = dr * sd;
    float A  = dot(q1, q1) - R1 * R1;
    float B  = dot(q0, q1) - R0 * R1;
    float C  = dot(q0, q0) - R0 * R0;
    if (abs(A) <= 1.0e-7) {
        // Along a generator, or the axis of a cylinder: 2 B t + C <= 0
        if (abs(B) < 1.0e-12) return C <= 0.0 ? vec2(-NO_HIT, NO_HIT) : none;
        float t = -0.5 * C / B;
        return B > 0.0 ? vec2(-NO_HIT, t) : vec2(t, NO_HIT);
    }
    float h = B * B - A * C;
    if (h < 0.0) return none;
    h = sqrt(h);
    // The roots without cancellation
    float q  = -(B + (B < 0.0 ? -h : h));
    float t1 = q / A;
    float t2 = q != 0.0 ? C / q : t1;
    vec2  r  = vec2(min(t1, t2), max(t1, t2));
    if (A > 0.0) {
        // Across the cone: within it between the roots, unless that is the other nappe, past the apex
        return R0 + R1 * 0.5 * (r.x + r.y) >= 0.0 ? r : none;
    }
    // Along the axis, within the half angle of the cone: within it beyond a root, on the side of the nappe (R > 0)
    return R1 > 0.0 ? vec2(r.y, NO_HIT) : vec2(-NO_HIT, r.x);
}

// The outward normal of the oblique cone at p (the gradient of the quadratic), also off it
vec3 oblique_cone_normal(vec3 p, vec3 pa, vec3 pb, float ra, float rb, vec3 w) {
    vec3  d  = pb - pa;
    float dw = dot(d, w);
    float dr = rb - ra;
    float s  = dot(p - pa, w) / dw;
    vec3  q  = p - pa - d * s;
    return q - w * ((dot(d, q) + (ra + dr * s) * dr) / dw);
}

// Where the ray enters the oblique cone of a side, cut to the segment (s, the interval between the planes of the rings,
// segment_interval, or where the chain ends the plane normal to w through the end of the side, the rung covers the
// rest): NO_HIT for none
float side_guess(vec3 ro, vec3 rd, vec3 pa, vec3 pb, float ra, float rb, vec3 w, vec2 s) {
    vec2 iv = oblique_cone_interval(ro, rd, pa, pb, ra, rb, w);
    if (!(iv.x <= iv.y)) return NO_HIT;
    if (in_frag.open_end.x == 1u) clip_to_plane(ro, rd, pa,  w, 0.0, s);
    if (in_frag.open_end.y == 1u) clip_to_plane(ro, rd, pb, -w, 0.0, s);
    float t0 = max(iv.x, s.x);
    return t0 <= min(iv.y, s.y) ? t0 : NO_HIT;
}

// Ray against a side as the mesh makes it, P(phi, v) = mix(pa + ra R0(phi), pb + rb R1(phi), v), R_i = cos(phi) X_i +
// sin(phi) Y_i (the frames of the rings, X = Y x Z), from where it enters the oblique cone at t_guess. The ray meets the
// line of P at phi where A = pa + ra R0 - ro, B = pb + rb R1 - ro and the ray are coplanar: f(phi) = (A x B).rd = 0, a
// trigonometric polynomial of degree 2. Then the point of that line on the ray, at v. Kept if it lies on the ray (to a
// fiftieth of r: at a grazing angle Newton's method converges slowly, and the guess is that close), within the segment
// (with a sliver of overlap, so no ray slips between the two at a joint) and nearer than t_hit.
void side_hit(vec3 ro, vec3 rd, vec3 pa, vec3 pb, float ra, float rb, vec3 X0, vec3 Y0, vec3 X1, vec3 Y1, vec3 w, float t_guess,
              inout float t_hit, inout vec3 n_hit, inout vec3 n_side) {
    vec3  d   = pb - pa;
    vec3  p   = ro + rd * t_guess;
    float v   = clamp(dot(p - pa, w) / dot(d, w), 0.0, 1.0);
    vec3  q   = p - mix(pa, pb, v);
    // (cos, sin) of the angle, as a unit vector
    vec2  e   = vec2(dot(q, mix(X0, X1, v)), dot(q, mix(Y0, Y1, v)));
    e = dot(e, e) > 0.0 ? normalize(e) : vec2(-1.0, 0.0);

    vec3  A0  = pa - ro;
    vec3  A1  = X0 * ra;
    vec3  A2  = Y0 * ra;
    vec3  B0  = pb - ro;
    vec3  B1  = X1 * rb;
    vec3  B2  = Y1 * rb;
    vec3  c0  = cross(B0, rd);
    vec3  c1  = cross(B1, rd);
    vec3  c2  = cross(B2, rd);
    // f = k0 + k1 cos + k2 sin + k11 cos^2 + k12 cos sin + k22 sin^2
    float k0  = dot(A0, c0);
    float k1  = dot(A1, c0) + dot(A0, c1);
    float k2  = dot(A2, c0) + dot(A0, c2);
    float k11 = dot(A1, c1);
    float k12 = dot(A1, c2) + dot(A2, c1);
    float k22 = dot(A2, c2);
    // Newton's steps, each along the tangent of the circle and back onto it (the angle moves by atan(step), which keeps
    // the convergence quadratic, without a trigonometric function)
    for (int i = 0; i < 2; ++i) {
        float f  = k0 + e.x * k1 + e.y * k2 + e.x * e.x * k11 + e.x * e.y * k12 + e.y * e.y * k22;
        float fp = e.x * k2 - e.y * k1 + 2.0 * e.x * e.y * (k22 - k11) + (e.x * e.x - e.y * e.y) * k12;
        if (fp != 0.0) e = normalize(e + vec2(e.y, -e.x) * (f / fp));
    }
    float cs = e.x;
    float sn = e.y;
    vec3  A  = A0 + A1 * cs + A2 * sn;
    vec3  E  = B0 + B1 * cs + B2 * sn - A;
    vec3  cE = cross(E, rd);
    float ee = dot(cE, cE);
    if (!(ee > 0.0)) return;
    v = -dot(cross(A, rd), cE) / ee;
    vec3  pt = A + E * v;
    float t  = dot(pt, rd);
    vec3  off = pt - rd * t;
    float tol = 0.02 * max(ra, rb);
    if (!(dot(off, off) <= tol * tol && v >= -1.0e-3 && v <= 1.0 + 1.0e-3 && t > 0.0 && t < t_hit)) return;

    // The normal, P_phi x P_v, outward (along the radius of the circles)
    vec3 R0 = X0 * cs + Y0 * sn;
    vec3 R1 = X1 * cs + Y1 * sn;
    vec3 n  = cross(mix((Y0 * cs - X0 * sn) * ra, (Y1 * cs - X1 * sn) * rb, v), E);
    if (dot(n, mix(R0 * ra, R1 * rb, v)) < 0.0) n = -n;
    t_hit  = t;
    n_hit  = n;
    n_side = n;
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

    // The faces the ray can enter (see above)
    vec3 ca = a1 - a0;
    vec3 cb = b1 - b0;
    float rs = r0 + r1;
    bool band    = dot(ca, ca) + dot(cb, cb) > 1.0e-8 * rs * rs;
    bool face_in = false;
    if (band) {
        vec3  y0 = in_frag.normal[0] * r0;
        vec3  y1 = in_frag.normal[1] * r1;
        vec3  ym = normalize(in_frag.normal[0] + in_frag.normal[1]);
        float sy = dot(rd, ym);
        float t0 = t;
        if (sy <  in_frag.radius.z) patch_hit(ro, rd, a0 + y0, a1 + y0, b0 + y1, b1 + y1, t, n, ns);
        bool  top = t < t0;
        float t1 = t;
        if (sy > -in_frag.radius.z) patch_hit(ro, rd, a0 - y0, a1 - y0, b0 - y1, b1 - y1, t, n, ns);
        if (t < t1) top = false;
        // Into the band through a face (against its outward normal, which leans to +-ym): the sides lie within the faces
        // (each line of a side is that of the face moved towards the core, and out past its end), a ray which enters
        // a face from outside has met neither before
        face_in = t < NO_HIT && dot(rd, n) * dot(n, top ? ym : -ym) < 0.0;
    }

    // A side, or the whole of a round tube; the other side: the one the guess meets first, then the other where it may
    // still be nearer. Where the guess is behind the nearest hit so far by more than the guess can be off (a tenth of
    // r), the side is behind it too.
    if (!face_in) {
        vec3  wt  = in_frag.tangent[0] + in_frag.tangent[1];
        vec3  X0  = cross(in_frag.normal[0], in_frag.tangent[0]);
        vec3  X1  = cross(in_frag.normal[1], in_frag.tangent[1]);
        vec2  seg = segment_interval(ro, rd);
        vec3  w0  = side_axis(b0 - a0, wt);
        vec3  w1  = band ? side_axis(b1 - a1, wt) : w0;
        float g0  = side_guess(ro, rd, a0, b0, r0, r1, w0, seg);
        float g1  = band ? side_guess(ro, rd, a1, b1, r0, r1, w1, seg) : NO_HIT;
        bool  near1 = g1 < g0;
        float slack = 0.1 * rs;
        if (!band) {
            // A round tube is the oblique cone (see above)
            if (g0 < t) {
                t  = g0;
                n  = oblique_cone_normal(ro + rd * g0, a0, b0, r0, r1, w0);
                ns = n;
            }
        } else {
            for (int k = 0; k < 2; ++k) {
                bool  one = (k == 0) == near1;
                float tg  = one ? g1 : g0;
                if (tg < NO_HIT && tg < t + slack) {
                    side_hit(ro, rd, one ? a1 : a0, one ? b1 : b0, r0, r1, X0, in_frag.normal[0], X1, in_frag.normal[1], one ? w1 : w0, tg, t, n, ns);
                }
            }
        }
    }

    // Where the chain ends, the rungs (a sphere for a round tube)
    if (in_frag.open_end.x == 1u) round_cone_hit(ro, rd, a0, a1, r0, r0, t, n, ns);
    if (in_frag.open_end.y == 1u) round_cone_hit(ro, rd, b0, b1, r1, r1, t, n, ns);

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
