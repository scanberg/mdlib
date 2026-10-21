#version 330 core

// Bonds (licorice, weak bonds, and the dashed partners which show delocalized bonds) as raycast impostors,
// drawn by vertex pulling without a geometry shader.
//
// Every primitive is one quad of 4 vertices: primitive p (quad_primitive, QUADS in md_gl.c), corner c = gl_VertexID & 3.
// The quads are drawn in blocks from a static index buffer (md_gl.c), two triangles (0, 1, 2), (2, 1, 3) per quad, each
// block once per instance of the molecule (INSTANCES): the atoms are taken into world space by its transform. A
// primitive is a bond, or for BOND_PARTNERS an entry of the partners.
//
// The quad bounds the image of the capsule around the segment (all points within the radius of it), which contains
// the solid cylinder and its dashes. It is a trapezoid in the image plane:
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
//
// BOND ORDERS (md_gl_mol_set_bond_orders). The plane of a bond is given by its surroundings: atoms bonded to either end,
// which lie in it.
//   - Double and triple bonds (covalent bonds, u_bond_orders != 0) are two or three thinner sticks in place of the
//     stick, side by side in the plane of the sp2 centers (bond_side): up to two surrounding atoms per end, each end
//     weighed by how planar it is in the frame shown. A triple bond, or a double bond without a plane (O=O, the C=O of
//     CO2, one between pyramidal ends), has them perpendicular to the view. The thinner sticks stay within u_spread
//     times the radius of the stick, and the bond keeps the spheres of the radius of the stick at its ends: with
//     u_spread 1 (licorice) the thinner sticks end within them, as within the joint the representation has at an atom.
//     The quad bounds the capsule of the axis which contains all of it.
//   - Delocalized bonds keep their stick and get a thin dashed partner (BOND_PARTNERS), a primitive of its own: an entry
//     is a bond and the atoms which give the side of the partner, within its group. So does an aromatic bond in no
//     circle (aromatic rings are circles, ring.vert), on the side of the aromatic bonds that continue from it. The ends
//     of the partner are inset to where the partners of the neighboring bonds at the same offset meet, as on an inner
//     polygon, and a little more. The quad bounds the partner alone.
//
// CHAIN GAPS (md_gl_mol_set_chain_gaps). A primitive is an entry of the gaps: (a, neighbour of a, b, neighbour of b). The
// gap runs from 5/6 of the way from the neighbour of a to a, where the cartoon of a chain ends (its uniform cubic
// B-spline at the clamped end of a chain), to the same point at b; an end which is its own neighbour is the atom itself.
//
// DASHED LINES (weak bonds, partners, gaps): dashes of u_dash_length with u_gap_length between them (Ångström), never a
// number of dashes. A line of length L gets the number of dashes n that comes closest to the lengths,
// n = round((L + gap) / (dash + gap)), at least 2, and bond.frag draws the first from pa and the last up to pb, the
// dashes and gaps stretched or shrunk alike by L / (n dash + (n - 1) gap) to fit: no dash is cut off at either end.

#ifndef ORTHO
#define ORTHO 0
#endif

// 1: weak bonds of the representation, with a weight per bond
#ifndef WEAK_BONDS
#define WEAK_BONDS 0
#endif

// 1: the dashed partners of delocalized covalent bonds (and of aromatic bonds in no circle)
#ifndef BOND_PARTNERS
#define BOND_PARTNERS 0
#endif

// 1: the gaps of the chains of the molecule, dashed lines over the breaks
#ifndef CHAIN_GAPS
#define CHAIN_GAPS 0
#endif

// Pixels the sides of the quad are moved out by
#define SUBPIXEL_MARGIN (1.0 / 32.0)

#define REF_NONE        0xFFFFFFFFu

// Kinds of bonds (md_gl.c, bond.frag). An elastic bond is drawn as a plain one, where it is drawn at all.
#define KIND_PLAIN      0u
#define KIND_DOUBLE     1u
#define KIND_TRIPLE     2u
#define KIND_ELASTIC    3u

// Relative to the radius of the stick (r) and to the envelope of the thinner sticks (u_spread * r)
#define DOUBLE_RADIUS   0.42    // * envelope
#define TRIPLE_RADIUS   0.28    // * envelope
#define PARTNER_RADIUS  0.33    // * r
#define PARTNER_OFFSET  1.9     // * envelope, axis to axis
#define PARTNER_GAP     0.08    // * length of the bond, added to the inset of either end

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
    float u_dash_length;
    float u_gap_length;
    uint  u_bond_orders;
    float u_spread;
    uint  u_gap_index_base;
    uint  u_bond_kinds;     // u_buf_bond_order is bound: the kind of every bond is known
    uint  u_elastic_bonds;  // Draw the elastic bonds (KIND_ELASTIC), which are left out otherwise
};

#pragma INSTANCES

uniform vec2 u_viewport_size;           // Pixels

uniform samplerBuffer  u_buf_atom_pos;      // RGB32F
uniform samplerBuffer  u_buf_atom_vel;      // RGB32F
uniform usamplerBuffer u_buf_atom_flags;    // R8UI
uniform samplerBuffer  u_buf_atom_color;    // RGBA8
#if CHAIN_GAPS
uniform usamplerBuffer u_buf_gap;           // RGBA32UI, (a, neighbour of a, b, neighbour of b)
uniform usamplerBuffer u_buf_gap_comp;      // R32UI, the component before the break, which the gap picks as
#else
uniform usamplerBuffer u_buf_bond;          // RG32UI, atom indices
#endif
#if WEAK_BONDS
uniform samplerBuffer  u_buf_bond_weight;   // R32F, scales the radius
#elif CHAIN_GAPS
#elif BOND_PARTNERS
uniform usamplerBuffer u_buf_bond_partner;  // RGBA32UI, (bond, atom bonded to its first atom, to its second atom, -)
#else
uniform usamplerBuffer u_buf_bond_order;    // R32UI, kind | index into u_buf_bond_ref << 3
uniform usamplerBuffer u_buf_bond_ref;      // RGBA32UI, two atoms bonded to the first atom, two to the second, or REF_NONE
#endif

out Fragment {
    flat vec3  pa;
    flat vec3  pb;
    flat float radius;
    flat vec4  color[2];
    flat vec3  view_vel[2];
    flat uint  atom_picking_idx[2];
    flat uint  bond_picking_idx;
    flat uint  flags[2];                // Of the atoms, a half each
#if BOND_PARTNERS
    flat vec2  seg_range;               // The part of the bond the partner spans, in [0, 1] from its first atom
#elif !WEAK_BONDS && !CHAIN_GAPS
    flat uint  kind;
    flat vec3  side;                    // Unit, perpendicular to the axis, in the plane of the bond
    flat vec2  sub;                     // (offset, radius) of the thinner sticks
#endif
#if WEAK_BONDS || BOND_PARTNERS || CHAIN_GAPS
    flat uint  dash_count;              // DASHED LINES
#endif
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

// The surroundings of a bond (pa, pb, view space): the atoms ref.x bonded to pa and ref.y bonded to pb, as vectors from
// their end (ca, cb). One farther than a bond can be is across the periodic boundary, and does not count.
void load_surroundings(uvec2 ref, vec3 pa, vec3 pb, out vec3 ca, out vec3 cb, out bool has_a, out bool has_b) {
    ca = vec3(0.0);
    cb = vec3(0.0);
    has_a = false;
    has_b = false;
    if (ref.x != REF_NONE) {
        ca = vec3(u_world_to_view * vec4(instance_point(texelFetch(u_buf_atom_pos, int(ref.x)).xyz), 1.0)) - pa;
        has_a = dot(ca, ca) < u_max_d2;
    }
    if (ref.y != REF_NONE) {
        cb = vec3(u_world_to_view * vec4(instance_point(texelFetch(u_buf_atom_pos, int(ref.y)).xyz), 1.0)) - pb;
        has_b = dot(cb, cb) < u_max_d2;
    }
}

// Unit direction perpendicular to the axis ax and to the view of the bond (pa, pb)
vec3 view_side(vec3 ax, vec3 pa, vec3 pb) {
#if ORTHO
    vec3 v = vec3(0.0, 0.0, 1.0);
#else
    vec3 v = -0.5 * (pa + pb);
#endif
    vec3 s = cross(ax, v);
    if (dot(s, s) <= 1.0e-12 * dot(v, v)) {
        s = cross(ax, abs(ax.x) < 0.9 ? vec3(1.0, 0.0, 0.0) : vec3(0.0, 1.0, 0.0));
    }
    return normalize(s);
}

// Unit direction in the plane of the bond (axis ax), perpendicular to the axis, for the side of a dashed partner. The
// normal of the plane is the sum of the normals the surroundings give, aligned. A surrounding atom within 6 degrees of
// the axis gives none; without a plane, perpendicular to the view.
vec3 plane_side(vec3 ax, vec3 pa, vec3 pb, vec3 ca, vec3 cb, bool has_a, bool has_b) {
    vec3  n  = vec3(0.0);
    float c2 = 0.0;
    if (has_a) {
        n = cross(ax, ca);
        c2 = dot(ca, ca);
    }
    if (has_b) {
        vec3 nb = cross(ax, cb);
        n += dot(n, nb) < 0.0 ? -nb : nb;
        c2 = max(c2, dot(cb, cb));
    }
    if (c2 > 0.0 && dot(n, n) > 0.01 * c2) {
        return normalize(cross(n, ax));
    }
    return view_side(ax, pa, pb);
}

// The weight of the side perpendicular to the view (bond_side): about that of one surrounding atom 8 degrees off the
// axis, sin^2
#define VIEW_SIDE_WEIGHT 0.02

// A surrounding atom (REF_NONE for none) of the end at p, as a vector c from it. One farther than a bond can be is
// across the periodic boundary, and does not count.
bool surrounding(uint ref, vec3 p, out vec3 c) {
    c = vec3(0.0);
    if (ref == REF_NONE) return false;
    c = vec3(u_world_to_view * vec4(instance_point(texelFetch(u_buf_atom_pos, int(ref)).xyz), 1.0)) - p;
    return dot(c, c) < u_max_d2;
}

// The line of c perpendicular to the axis, as its doubled angle in the basis (u, v) of the plane perpendicular to the
// axis, of length sin^2 of the angle between c and the axis: lines of opposite directions are the same
vec2 doubled_angle(vec3 c, vec3 u, vec3 v) {
    float x = dot(c, u);
    float y = dot(c, v);
    return vec2(x * x - y * y, 2.0 * x * y) / dot(c, c);
}

// The lines the surroundings (ref, up to two) of the end at p give, as doubled angles (doubled_angle). Two count as
// much as the end is planar: fully while the bond is within 10 degrees of their plane, not at all past 30. An end that
// is pyramidal in the frame shown has no plane of its own (an sp3 center, whatever the order of the bond says).
vec2 end_lines(vec3 ax, vec3 p, uvec2 ref, vec3 u, vec3 v) {
    vec3 c0, c1;
    bool has0 = surrounding(ref.x, p, c0);
    bool has1 = surrounding(ref.y, p, c1);
    vec2 s = vec2(0.0);
    if (has0) s += doubled_angle(c0, u, v);
    if (has1) s += doubled_angle(c1, u, v);
    if (has0 && has1) {
        // sin of the angle between the bond and the plane of the two (0 for two in a line, whose plane any is)
        vec3  m = cross(c0, c1);
        float sin_out = abs(dot(ax, m)) * inversesqrt(max(dot(m, m), 1.0e-20));
        s *= 1.0 - smoothstep(0.17, 0.5, sin_out);
    }
    return s;
}

// Unit direction of the thinner sticks of a double or triple bond (axis ax, ends pa and pb), perpendicular to the axis:
// in the plane of the sp2 centers. Its sign is of no consequence, the sticks lie symmetrically about the axis. ref: two
// surrounding atoms of either end (md_gl.c, bo_plane_refs), none for a triple bond.
// Each surrounding atom gives a line, the plane of the bond and the atom, and the lines are averaged as such: as
// doubled angles about the axis, weighted by sin^2 of the angle to the axis (end_lines). So cis and trans substituents
// agree, an end twisted against the other moves the side smoothly instead of flipping it, an atom close to the axis
// fades out, and a pyramidal end gives nothing. The side perpendicular to the view is a weak line of its own: it takes
// over where the surroundings give none (a triple bond, O=O, the C=O of CO2), and the plane fades into it rather than
// switching where they give little.
vec3 bond_side(vec3 ax, vec3 pa, vec3 pb, uvec4 ref) {
    vec3 u = view_side(ax, pa, pb);
    vec3 v = cross(ax, u);
    vec2 s = vec2(VIEW_SIDE_WEIGHT, 0.0) + end_lines(ax, pa, ref.xy, u, v) + end_lines(ax, pb, ref.zw, u, v);
    // Half the angle of s: (cos, sin) of it is along (|s| + s.x, s.y) and along (s.y, |s| - s.x), the larger of the two
    float l = length(s);
    vec2  h = s.x >= 0.0 ? vec2(l + s.x, s.y) : vec2(s.y, l - s.x);
    float hh = dot(h, h);
    if (!(hh > 0.0)) return u;
    h *= inversesqrt(hh);
    return h.x * u + h.y * v;
}

// 1 / tan(angle / 2) for the angle at an end between the bond (u, unit) and its surrounding atom (c): where the
// parallels to the two bonds at the same distance on the inside meet, along the bond, per unit of distance. 1 / tan(60)
// without one.
float inset_factor(vec3 u, vec3 c, bool has) {
    if (!has) return 0.57735;
    float cos_angle = clamp(dot(u, c) * inversesqrt(dot(c, c)), -0.999, 0.999);
    return sqrt((1.0 + cos_angle) / (1.0 - cos_angle));
}

void main() {
    uint prim   = quad_primitive();
    uint corner = uint(gl_VertexID) & 3u;
    if (prim == QUAD_NONE) {
        cull();
        return;
    }

#if BOND_PARTNERS
    uvec4 partner = texelFetch(u_buf_bond_partner, int(prim));
    uint bond = partner.x;
#else
    uint bond = prim;
#endif

#if !WEAK_BONDS && !BOND_PARTNERS && !CHAIN_GAPS
    // ELASTIC BONDS (MD_BOND_KIND_ELASTIC): springs of the force field, not bonds of the molecule, drawn only where the
    // representation shows them
    uint order = u_bond_kinds != 0u ? texelFetch(u_buf_bond_order, int(bond)).x : KIND_PLAIN;
    if ((order & 7u) == KIND_ELASTIC && u_elastic_bonds == 0u) {
        cull();
        return;
    }
#endif

#if CHAIN_GAPS
    uvec4 gap  = texelFetch(u_buf_gap, int(prim));
    uvec2 atom = gap.xz;
#else
    uvec2 atom = texelFetch(u_buf_bond, int(bond)).xy;
#endif

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

    vec3 wa = texelFetch(u_buf_atom_pos, int(atom.x)).xyz;
    vec3 wb = texelFetch(u_buf_atom_pos, int(atom.y)).xyz;
    vec3 va = texelFetch(u_buf_atom_vel, int(atom.x)).xyz;
    vec3 vb = texelFetch(u_buf_atom_vel, int(atom.y)).xyz;
#if CHAIN_GAPS
    // Where the cartoon ends: 5/6 of the way from the neighbour (the atom itself for a link)
    wa = (5.0 * wa + texelFetch(u_buf_atom_pos, int(gap.y)).xyz) * (1.0 / 6.0);
    wb = (5.0 * wb + texelFetch(u_buf_atom_pos, int(gap.w)).xyz) * (1.0 / 6.0);
    va = (5.0 * va + texelFetch(u_buf_atom_vel, int(gap.y)).xyz) * (1.0 / 6.0);
    vb = (5.0 * vb + texelFetch(u_buf_atom_vel, int(gap.w)).xyz) * (1.0 / 6.0);
#endif
    vec3 pa = vec3(u_world_to_view * vec4(instance_point(wa), 1.0));
    vec3 pb = vec3(u_world_to_view * vec4(instance_point(wb), 1.0));

    vec3  d  = pb - pa;
    float d2 = dot(d, d);
    if (!(r > 0.0) || d2 <= 0.0 || u_max_d2 <= d2) {
        cull();
        return;
    }

    out_frag.flags[0] = flags_a;
    out_frag.flags[1] = flags_b;
    out_frag.color[0] = color_a;
    out_frag.color[1] = color_b;
    out_frag.view_vel[0] = mat3(u_world_to_view) * instance_velocity(wa, va);
    out_frag.view_vel[1] = mat3(u_world_to_view) * instance_velocity(wb, vb);
    out_frag.atom_picking_idx[0] = u_atom_base_index + atom.x;
    out_frag.atom_picking_idx[1] = u_atom_base_index + atom.y;
#if CHAIN_GAPS
    out_frag.bond_picking_idx    = u_gap_index_base + texelFetch(u_buf_gap_comp, int(prim)).x;
#else
    out_frag.bond_picking_idx    = u_bond_base_index + bond;
#endif

#if BOND_PARTNERS
    {
        // The partner, on the side of the surroundings of this entry, in their plane
        float len = sqrt(d2);
        vec3  ax  = d / len;
        vec3  ca, cb;
        bool  has_a, has_b;
        load_surroundings(partner.yz, pa, pb, ca, cb, has_a, has_b);
        vec3 side = plane_side(ax, pa, pb, ca, cb, has_a, has_b);
        if ((has_a && dot(side, ca) < 0.0) || (!has_a && has_b && dot(side, cb) < 0.0)) {
            side = -side;
        }
        float fa = inset_factor( ax, ca, has_a);
        float fb = inset_factor(-ax, cb, has_b);
        // Not past a fraction of the inradius of the ring, len / 2 * tan(angle / 2) for a regular one
        float offset = min(PARTNER_OFFSET * u_spread * r, 0.225 * len / max(fa, fb));
        float ia = offset * fa + PARTNER_GAP * len;
        float ib = offset * fb + PARTNER_GAP * len;
        float k  = min(1.0, 0.8 * len / (ia + ib));
        ia *= k;
        ib *= k;
        out_frag.seg_range = vec2(ia / len, 1.0 - ib / len);
        // From here on the primitive is the partner
        pa = pa + side * offset + ax * ia;
        pb = pb + side * offset - ax * ib;
        d  = pb - pa;
        d2 = dot(d, d);
        r  = PARTNER_RADIUS * r;
    }
    out_frag.pa = pa;
    out_frag.pb = pb;
    out_frag.radius = r;
#elif WEAK_BONDS
    out_frag.pa = pa;
    out_frag.pb = pb;
    out_frag.radius = r;
#elif CHAIN_GAPS
    out_frag.pa = pa;
    out_frag.pb = pb;
    out_frag.radius = r;
#else
    out_frag.pa = pa;
    out_frag.pb = pb;
    out_frag.radius = r;

    // Without the orders shown every bond is a single stick, and so is an elastic one where it is shown
    uint kind = u_bond_orders != 0u ? (order & 7u) : KIND_PLAIN;
    vec3 side = vec3(0.0);
    vec2 sub  = vec2(0.0);
    if (kind == KIND_DOUBLE || kind == KIND_TRIPLE) {
        vec3 ax = d * inversesqrt(d2);
        side = bond_side(ax, pa, pb, texelFetch(u_buf_bond_ref, int(order >> 3u)));
        float envelope = u_spread * r;
        float rs = (kind == KIND_DOUBLE ? DOUBLE_RADIUS : TRIPLE_RADIUS) * envelope;
        sub = vec2(envelope - rs, rs);
        // The quad bounds the thinner sticks and the spheres of radius r at the ends (bond.frag)
        r = max(envelope, r);
    } else {
        kind = KIND_PLAIN;
    }
    out_frag.kind = kind;
    out_frag.side = side;
    out_frag.sub  = sub;
#endif

#if WEAK_BONDS || BOND_PARTNERS || CHAIN_GAPS
    {
        // DASHED LINES: the number of dashes which comes closest to their length and that of the gaps, over the line
        // as drawn (a partner inset from its bond)
        float dash = max(u_dash_length, 1.0e-4);
        float gap  = max(u_gap_length, 0.0);
        float n = floor((sqrt(d2) + gap) / (dash + gap) + 0.5);
        out_frag.dash_count = uint(clamp(n, 2.0, 4096.0));
    }
#endif

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
