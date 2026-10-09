#version 330 core

// Bonds (licorice, weak bonds) as raycast impostors, drawn by vertex pulling without a geometry shader.
//
// Every bond is one quad of 4 vertices: bond b, corner c = gl_VertexID & 3 with b = u_bond_offset +
// gl_InstanceID * u_quads_per_instance + (gl_VertexID >> 2). The quads are drawn from a static index buffer (md_gl.c),
// two triangles (0, 1, 2), (2, 1, 3) per quad.
//
// The quad bounds the image of the capsule around the bond (all points within the radius of the segment), which
// contains the solid cylinder and its dashes. It is a trapezoid in the image plane:
//   - The long sides are the two planes through the eye that are tangent to the capsule. They are parallel to the
//     axis and touch the infinite cylinder of the bond, so they are the exact outline of its body.
//   - The short sides are perpendicular to the bisector of the long sides and tangent to the outline of the end
//     spheres. The capsule lies in the wedge of the long sides, which opens away from its apex (the vanishing point of
//     the axis), so the short sides cut the wedge before the apex and the trapezoid is convex.
// The capsule is the convex hull of its end spheres, and a central projection keeps convex hulls for points in front
// of the eye, so the trapezoid bounds its image. Where the trapezoid is larger than the rectangle of the same
// orientation bounding both spheres (the eye close to the axis, where the long sides diverge), the rectangle is used
// instead, as it is when there are no tangent planes (the eye within the radius of the axis).
// In orthographic projection the image is a stadium and the quad is its bounding rectangle, which is exact.
//
// Only the part of the capsule behind the near plane is visible. A capsule that crosses it is bounded by the box of
// that part instead of the end spheres, since it may extend behind the eye, where its image is unbounded.
//
// The quad is placed at the depth of the point of the capsule closest to the image plane, or at the near plane if that
// is in front of it. It is in front of every visible point of the capsule, so the raycast starts at the quad, and its
// depth is a lower bound of the depth the fragment shader writes (layout (depth_greater) in bond.frag).

#ifndef ORTHO
#define ORTHO 0
#endif

// 0: Covalent bonds of the molecule (licorice), 1: weak bonds of the representation, with a weight per bond
#ifndef WEAK_BONDS
#define WEAK_BONDS 0
#endif

// Pixels the sides of the quad are moved out by
#define SUBPIXEL_MARGIN (1.0 / 32.0)

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

uniform uint u_bond_offset;
uniform uint u_quads_per_instance;
uniform vec2 u_viewport_size;           // Pixels

uniform samplerBuffer  u_buf_atom_pos;      // RGB32F
uniform samplerBuffer  u_buf_atom_vel;      // RGB32F
uniform usamplerBuffer u_buf_atom_flags;    // R8UI
uniform samplerBuffer  u_buf_atom_color;    // RGBA8
uniform usamplerBuffer u_buf_bond;          // RG32UI, atom indices
#if WEAK_BONDS
uniform samplerBuffer  u_buf_bond_weight;   // R32F, scales the radius
#endif

out Fragment {
    flat vec3  pa;
    flat vec3  pb;
    flat float radius;
    flat vec4  color[2];
    flat vec3  view_vel[2];
    flat uint  atom_picking_idx[2];
    flat uint  bond_picking_idx;
    smooth vec3 view_pos;               // On the quad, the ray of the fragment starts here
} out_frag;

// Support of the image of the sphere (c, r) along the unit direction n of the image plane z = -1: the lines n.x = s
// which are tangent to its outline, returned as (s_min, s_max). The plane through the eye and such a line has the
// normal (n, s), and it is tangent to the sphere where (n.c.xy + s c.z)^2 = r^2 (1 + s^2).
// Requires the sphere to be in front of the eye (c.z + r < 0).
vec2 sphere_support(vec3 c, float r, vec2 n) {
    float m = dot(n, c.xy);
    float k = c.z * c.z - r * r;
    float h = r * sqrt(max(k + m * m, 0.0));
    float s = -c.z * m;
    return vec2(s - h, s + h) / k;
}

// Support of the image of the part of the box [lo, hi] (view space) behind the near plane, along the unit direction n
// of the image plane z = -1, as (s_min, s_max). n.xy over the box is within the center +- the extent, and divided by
// the depth in [near, -lo.z] it is extremal at either end.
vec2 box_support(vec3 lo, vec3 hi, float near, vec2 n) {
    float m  = dot(n, 0.5 * (lo.xy + hi.xy));
    float e  = dot(abs(n), 0.5 * (hi.xy - lo.xy));
    float d0 = near;
    float d1 = max(-lo.z, near);
    return vec2(min((m - e) / d0, (m - e) / d1), max((m + e) / d0, (m + e) / d1));
}

void cull() {
    // All corners of the quad coincide: no area, no fragments
    gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
}

void main() {
    uint bond   = u_bond_offset + uint(gl_InstanceID) * u_quads_per_instance + (uint(gl_VertexID) >> 2u);
    uint corner = uint(gl_VertexID) & 3u;

    uvec2 atom = texelFetch(u_buf_bond, int(bond)).xy;

    uint flags_a = texelFetch(u_buf_atom_flags, int(atom.x)).x;
    uint flags_b = texelFetch(u_buf_atom_flags, int(atom.y)).x;
    vec4 color_a = texelFetch(u_buf_atom_color, int(atom.x));
    vec4 color_b = texelFetch(u_buf_atom_color, int(atom.y));

    if ((flags_a & u_atom_mask) != u_atom_mask || (flags_b & u_atom_mask) != u_atom_mask || color_a.a == 0.0 || color_b.a == 0.0) {
        cull();
        return;
    }

    float r = u_radius;
#if WEAK_BONDS
    r *= texelFetch(u_buf_bond_weight, int(bond)).x;
#endif

    vec3 pa = vec3(u_world_to_view * vec4(texelFetch(u_buf_atom_pos, int(atom.x)).xyz, 1.0));
    vec3 pb = vec3(u_world_to_view * vec4(texelFetch(u_buf_atom_pos, int(atom.y)).xyz, 1.0));

    vec3  d  = pb - pa;
    float d2 = dot(d, d);
    if (!(r > 0.0) || d2 <= 0.0 || u_max_d2 <= d2) {
        cull();
        return;
    }

    out_frag.pa = pa;
    out_frag.pb = pb;
    out_frag.radius = r;
    out_frag.color[0] = color_a;
    out_frag.color[1] = color_b;
    out_frag.view_vel[0] = mat3(u_world_to_view) * texelFetch(u_buf_atom_vel, int(atom.x)).xyz;
    out_frag.view_vel[1] = mat3(u_world_to_view) * texelFetch(u_buf_atom_vel, int(atom.y)).xyz;
    out_frag.atom_picking_idx[0] = u_atom_base_index + atom.x;
    out_frag.atom_picking_idx[1] = u_atom_base_index + atom.y;
    out_frag.bond_picking_idx    = u_bond_base_index + bond;

    // Distance from the eye to the near plane (z = -near)
#if ORTHO
    float near = (u_view_to_clip[3][2] + 1.0) / u_view_to_clip[2][2];
#else
    float near = u_view_to_clip[3][2] / (u_view_to_clip[2][2] - 1.0);
#endif

    // Extent of the capsule along z
    float z_front = max(pa.z, pb.z) + r;
    float z_back  = min(pa.z, pb.z) - r;
    if (z_back >= -near) {
        // Entirely in front of the near plane
        cull();
        return;
    }
    // The quad is in front of the part of the capsule which is behind the near plane: the part in front is clipped
    float z_quad = min(z_front, -near);

    // The sides of the quad are tangent to the outline. They are moved out by a fraction of a pixel (in units of the
    // image plane z = -1, or of the view plane in orthographic projection), so that neither rounding nor the snapping
    // of the vertices to the subpixel grid of the rasterizer loses fragments whose centers are just inside of it.
    vec2 px = 2.0 / (u_viewport_size * vec2(u_view_to_clip[0][0], u_view_to_clip[1][1]));
    float margin = u_viewport_size.x > 0.0 ? SUBPIXEL_MARGIN * max(abs(px.x), abs(px.y)) : 0.0;

    // Corner: (s, t) in the basis (g, gp) of the image plane
    bool hi = (corner & 1u) != 0u;
    bool up = (corner & 2u) != 0u;
    vec2 g, gp;
    float s, t;

#if ORTHO
    // The image is the stadium of the segment (pa.xy, pb.xy) and the radius
    vec2 a2 = pa.xy;
    vec2 b2 = pb.xy;
    vec2 ab = b2 - a2;
    float len = length(ab);
    g  = len > 1.0e-6 * r ? ab / len : vec2(1.0, 0.0);
    gp = vec2(-g.y, g.x);
    float sa = dot(a2, g);
    float sb = dot(b2, g);
    float ta = dot(a2, gp);
    float tb = dot(b2, gp);
    s = hi ? max(sa, sb) + r + margin : min(sa, sb) - r - margin;
    t = up ? max(ta, tb) + r + margin : min(ta, tb) - r - margin;
    vec3 view_pos = vec3(s * g + t * gp, z_quad);
#else
    // A capsule that crosses the near plane may extend behind the eye, where its image is unbounded. The part behind
    // the near plane lies within the capsule of the part of the axis behind z = -near + r (other points of the axis are
    // farther than r from it), and that within its box. The image of the box clipped at the near plane is bounded.
    bool crossing = z_front > -near;
    vec3 box_lo, box_hi;
    if (crossing) {
        float zc = -near + r;
        vec3 ca = pa;
        vec3 cb = pb;
        if (ca.z > zc) ca = mix(ca, cb, (ca.z - zc) / (ca.z - cb.z));
        if (cb.z > zc) cb = mix(cb, ca, (cb.z - zc) / (cb.z - ca.z));
        box_lo = min(ca, cb) - r;
        box_hi = max(ca, cb) + r;
    }

    vec3 ax = d * inversesqrt(d2);
    vec3 w  = pa - ax * dot(pa, ax);    // From the eye to the closest point of the axis
    float h2 = dot(w, w);

    bool trapezoid = false;
    vec2 n_lo, n_up;    // Unit inward normals of the long sides, n.x >= c, which bound t from below and above
    float c_lo, c_up;

    if (h2 > r * r * 1.0001) {
        // Planes through the eye at distance r from the axis and parallel to it: N.X = 0 with
        // N = (r e1 +- sqrt(h^2 - r^2) e2) / h. N.X = r on the axis, so the capsule is in N.X >= 0, which is N.xy.x >= N.z
        // in the image (X = (x, -1) * depth) for the points in front of the eye.
        float h  = sqrt(h2);
        vec3  e1 = w / h;
        vec3  e2 = cross(ax, e1);
        vec3  N0 = (r * e1 + sqrt(h2 - r * r) * e2) / h;
        vec3  N1 = (r * e1 - sqrt(h2 - r * r) * e2) / h;
        float l0 = length(N0.xy);
        float l1 = length(N1.xy);
        if (l0 > 1.0e-4 && l1 > 1.0e-4) {
            vec2 n0 = N0.xy / l0;
            vec2 n1 = N1.xy / l1;
            // The wedge n0.x >= c0, n1.x >= c1 opens along n0 + n1, away from its apex
            vec2 bis = n0 + n1;
            g  = dot(bis, bis) > 1.0e-8 ? normalize(bis) : vec2(-n0.y, n0.x);
            gp = vec2(-g.y, g.x);
            float g0 = dot(n0, gp);
            float g1 = dot(n1, gp);
            // One side bounds t from below (n.gp > 0), the other from above
            if (abs(g0) > 1.0e-4 && abs(g1) > 1.0e-4 && g0 * g1 < 0.0) {
                trapezoid = true;
                n_lo = g0 > 0.0 ? n0 : n1;
                n_up = g0 > 0.0 ? n1 : n0;
                c_lo = (g0 > 0.0 ? N0.z / l0 : N1.z / l1) - margin;
                c_up = (g0 > 0.0 ? N1.z / l1 : N0.z / l0) - margin;
            }
        }
    }

    if (!trapezoid) {
        // Rectangle along the image of the axis, or of the image plane for the box
        vec2 ab = crossing ? vec2(0.0) : pb.xy / -pb.z - pa.xy / -pa.z;
        float len2 = dot(ab, ab);
        g  = len2 > 1.0e-12 ? ab * inversesqrt(len2) : vec2(1.0, 0.0);
        gp = vec2(-g.y, g.x);
    }

    // Short sides tangent to the end spheres, or to the box
    vec2 sg, st;
    if (crossing) {
        sg = box_support(box_lo, box_hi, near, g);
        st = box_support(box_lo, box_hi, near, gp);
    } else {
        vec2 sa = sphere_support(pa, r, g);
        vec2 sb = sphere_support(pb, r, g);
        vec2 ta = sphere_support(pa, r, gp);
        vec2 tb = sphere_support(pb, r, gp);
        sg = vec2(min(sa.x, sb.x), max(sa.y, sb.y));
        st = vec2(min(ta.x, tb.x), max(ta.y, tb.y));
    }
    float s_lo = sg.x - margin;
    float s_hi = sg.y + margin;
    // Long sides of the rectangle
    float t_lo = st.x - margin;
    float t_hi = st.y + margin;

    if (trapezoid) {
        // The long sides at the short sides: n.(s g + t gp) = c
        float lo0 = (c_lo - s_lo * dot(n_lo, g)) / dot(n_lo, gp);
        float lo1 = (c_lo - s_hi * dot(n_lo, g)) / dot(n_lo, gp);
        float up0 = (c_up - s_lo * dot(n_up, g)) / dot(n_up, gp);
        float up1 = (c_up - s_hi * dot(n_up, g)) / dot(n_up, gp);
        // The width grows along g. The box (and rounding) may reach past the apex, which then becomes the narrow end.
        float w0 = up0 - lo0;
        float w1 = up1 - lo1;
        if (w0 < 0.0 && w1 > 0.0) {
            float f = -w0 / (w1 - w0);
            s_lo = mix(s_lo, s_hi, f);
            lo0  = mix(lo0, lo1, f);
            up0  = lo0;
            w0   = 0.0;
        }
        if (w0 >= 0.0 && (w0 + w1) * (s_hi - s_lo) < 2.0 * (t_hi - t_lo) * (sg.y - sg.x + 2.0 * margin)) {
            t = up ? (hi ? up1 : up0) : (hi ? lo1 : lo0);
        } else {
            trapezoid = false;
            s_lo = sg.x - margin;
        }
    }
    if (!trapezoid) {
        t = up ? t_hi : t_lo;
    }
    s = hi ? s_hi : s_lo;

    vec3 view_pos = vec3(s * g + t * gp, -1.0) * -z_quad;
#endif

    vec4 clip = u_view_to_clip * vec4(view_pos, 1.0);
    // Keep the quad behind the near plane, the depth does not change its image
    clip.z = max(clip.z, -clip.w * (1.0 - 1.0e-6));

    out_frag.view_pos = view_pos;
    gl_Position = clip;
}
