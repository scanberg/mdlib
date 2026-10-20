#version 330 core

// Aromatic rings (md_gl_mol_set_bond_orders) as raycast tori: a circle inside each ring, drawn with the bonds of the
// licorice and ball and stick representations. Vertex pulling as bond.vert, from the same index buffer: primitive
// p (quad_primitive, QUADS in md_gl.c) is one quad of 4 vertices, a ring, of an instance of the molecule (INSTANCES).
//
// The circle follows the atoms of the ring in every frame. In view space:
//   - center: the centroid of the atoms
//   - normal: Newell's, the sum of the cross products of consecutive atoms about the center. Robust for the puckering
//     of a ring in motion, and along the order of the atoms: they go counter clockwise about it.
//   - major radius: the apothem, the smallest distance in the plane from the center to an edge, less the offset the
//     dashed partners have from their bond (bond.vert), at most a fraction of the apothem. The circle sits where the
//     dashes of a delocalized group do, and keeps clear of the sticks in an irregular ring (a pyrrole).
//   - tube radius: a fraction of the radius of the sticks, as the dashes
//   - dashes: none (a solid circle) or the dashes over each bond of a dashed circle, the circle of a ring of a fused
//     aromatic system (md_gl.c). DASHED_PER_BOND, or one where those would be shorter than the tube is thick (thick
//     sticks), as in a ring of the same size with bonds of RING_IDEAL_BOND and the same sticks: the number follows from
//     the size of the ring and the sticks alone, and does not change from frame to frame as the atoms move.
// A ring is culled when one of its atoms is hidden (alpha 0, u_atom_mask), as a bond is when either atom is: so it is in
// the selection and highlight passes exactly when all of its atoms are. Also when an edge is longer than a bond can be
// (across the periodic boundary), or when the ring has no plane or no room inside for the circle.
//
// The quad bounds the image of the box of the torus, half extents R + r in its plane and r along its normal, and is
// placed at the depth of the point of the box closest to the image plane (or at the near plane, if that is in front of
// it). It is in front of every visible point of the torus, so the raycast starts at the quad, and its depth is a lower
// bound of the depth ring.frag writes (layout (depth_greater)).

#ifndef ORTHO
#define ORTHO 0
#endif

// Pixels the sides of the quad are moved out by
#define SUBPIXEL_MARGIN (1.0 / 32.0)

// Entries of the rings (md_gl.c, GL_BUFFER_AROMATIC_RING): two RGBA32UI texels, the index of the ring | dashed <<
// RING_DASHED_SHIFT | its atom count << RING_COUNT_SHIFT, then up to MAX_RING_ATOMS atoms in order around it
#define MAX_RING_ATOMS      7
#define RING_DASHED_SHIFT   28u
#define RING_COUNT_SHIFT    29u

// Relative to the radius of the sticks (r), as the dashed partners of bond.vert
#define RING_TUBE_RADIUS    0.33    // * r
#define RING_OFFSET         1.9     // * u_spread * r, from the axes of the bonds to the axis of the tube
#define RING_MAX_INSET      0.45    // * apothem, the most the offset may take...
#define RING_MIN_GAP        0.5     // ...unless the tube would then come closer to a stick than this * its radius

// Dashed circles (ring.frag): dashes over each bond, the gaps centred on the atoms
#define DASHED_PER_BOND     2u
#define RING_DASH_FILL      0.5     // Of each dash and its gap, the dash (as ring.frag)
#define RING_IDEAL_BOND     1.39    // Å, of the ring the number of dashes is decided on
#define RING_MIN_DASH       1.0     // * the diameter of the tube, the shortest dash (otherwise one dash per bond)

#define PI 3.14159265358979

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

    float u_radius;             // Of the sticks
    float u_max_d2;             // Longest bond, squared
    int   u_mode;
    float u_sharpness;
    uint  u_uniform_color;
    float u_spread;
    uint  u_ring_base_index;    // Picking: a ring writes u_ring_base_index + its index
};

#pragma INSTANCES

uniform vec2 u_viewport_size;           // Pixels

uniform samplerBuffer  u_buf_atom_pos;      // RGB32F
uniform samplerBuffer  u_buf_atom_vel;      // RGB32F
uniform usamplerBuffer u_buf_atom_flags;    // R8UI
uniform samplerBuffer  u_buf_atom_color;    // RGBA8
uniform usamplerBuffer u_buf_ring;          // RGBA32UI, two per ring

out Fragment {
    flat vec3  center;          // View space
    flat vec3  normal;          // Unit, the atoms go counter clockwise about it
    flat vec3  tangent;         // Unit, in the plane towards the first atom: the x axis of the frame of the torus
    flat vec2  radii;           // Major, tube
    flat vec4  angle[2];        // Of atom k (in angle[k / 4][k % 4]) about the normal from the tangent, increasing from 0
    flat uvec4 color[2];        // Packed RGBA8 of atom k (in color[k / 4][k % 4])
    flat uvec2 count_picking;   // Atoms, picking index
    flat vec3  view_vel;        // Mean velocity of the atoms
    flat uint  flags;           // Of all its atoms: those they have in common
    flat uint  dashes;          // Over each bond, 0 for a solid circle
    smooth vec3 view_pos;       // On the quad, the ray of the fragment starts here
} out_frag;

// Support of the image of the part of the box [lo, hi] (view space) behind the near plane, along the unit direction n
// of the image plane z = -1, as (s_min, s_max). As in bond.vert.
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

uint pack_color(vec4 c) {
    uvec4 u = uvec4(round(clamp(c, 0.0, 1.0) * 255.0));
    return u.x | (u.y << 8u) | (u.z << 16u) | (u.w << 24u);
}

void main() {
    uint prim   = quad_primitive();
    uint corner = uint(gl_VertexID) & 3u;
    if (prim == QUAD_NONE) {
        cull();
        return;
    }

    uvec4 e0 = texelFetch(u_buf_ring, int(2u * prim + 0u));
    uvec4 e1 = texelFetch(u_buf_ring, int(2u * prim + 1u));
    uint ring_idx = e0.x & ((1u << RING_DASHED_SHIFT) - 1u);
    bool dashed   = ((e0.x >> RING_DASHED_SHIFT) & 1u) != 0u;
    int  n        = int(e0.x >> RING_COUNT_SHIFT);
    uint atom[MAX_RING_ATOMS] = uint[MAX_RING_ATOMS](e0.y, e0.z, e0.w, e1.x, e1.y, e1.z, e1.w);

    if (n < 3 || n > MAX_RING_ATOMS) {
        cull();
        return;
    }

    vec3 p[MAX_RING_ATOMS];
    uint col[MAX_RING_ATOMS];
    vec3 c   = vec3(0.0);
    vec3 vel = vec3(0.0);
    vec3 local_center = vec3(0.0);  // The centroid in the molecule, for the velocity
    uint common_flags = 0xFFFFFFFFu;
    for (int k = 0; k < MAX_RING_ATOMS; ++k) {
        col[k] = 0u;
        if (k >= n) continue;
        uint flags  = texelFetch(u_buf_atom_flags, int(atom[k])).x;
        vec4 color  = texelFetch(u_buf_atom_color, int(atom[k]));
        if ((flags & u_atom_mask) != u_atom_mask || color.a == 0.0) {
            cull();
            return;
        }
        common_flags &= flags;
        vec3 w = texelFetch(u_buf_atom_pos, int(atom[k])).xyz;
        p[k]   = vec3(u_world_to_view * vec4(instance_point(w), 1.0));
        local_center += w;
        col[k] = pack_color(color);
        vel   += texelFetch(u_buf_atom_vel, int(atom[k])).xyz;
        c     += p[k];
    }
    c /= float(n);

    // The edges are bonds: one longer than a bond can be is across the periodic boundary
    vec3 N = vec3(0.0);
    for (int k = 0; k < n; ++k) {
        vec3 a = p[k];
        vec3 b = p[(k + 1) % n];
        vec3 e = b - a;
        if (u_max_d2 <= dot(e, e)) {
            cull();
            return;
        }
        N += cross(a - c, b - c);
    }
    float N2 = dot(N, N);
    float size2 = dot(p[0] - c, p[0] - c);
    if (!(N2 > 1.0e-8 * size2 * size2)) {
        cull();
        return;
    }
    vec3 nrm = N * inversesqrt(N2);

    // Apothem: the smallest distance in the plane from the center to the line of an edge
    float apothem = 3.4e38;
    for (int k = 0; k < n; ++k) {
        vec3 a = p[k] - c;
        vec3 b = p[(k + 1) % n] - c;
        a -= nrm * dot(a, nrm);
        b -= nrm * dot(b, nrm);
        vec3 d = b - a;
        float dd = dot(d, d);
        if (dd > 0.0) {
            apothem = min(apothem, abs(dot(cross(a, d), nrm)) * inversesqrt(dd));
        }
    }

    // Thick sticks in a small ring: the circle keeps clear of them, or it is culled for want of room below
    float r_stick = u_radius;
    float r       = RING_TUBE_RADIUS * r_stick;
    float offset  = max(min(RING_OFFSET * u_spread * r_stick, RING_MAX_INSET * apothem), r_stick + (1.0 + RING_MIN_GAP) * r);
    float R = apothem - offset;
    if (!(r > 0.0) || !(R > 2.0 * r)) {
        cull();
        return;
    }

    // Dashes over each bond, from the ring of the same size with ideal bonds (in view space: scaled as the sticks are)
    uint dashes = 0u;
    if (dashed) {
        float scale   = (length(u_world_to_view[0].xyz) + length(u_world_to_view[1].xyz) + length(u_world_to_view[2].xyz)) / 3.0;
        float ideal_a = RING_IDEAL_BOND * scale / (2.0 * tan(PI / float(n)));
        float ideal_o = max(min(RING_OFFSET * u_spread * r_stick, RING_MAX_INSET * ideal_a), r_stick + (1.0 + RING_MIN_GAP) * r);
        float arc     = 2.0 * PI * (ideal_a - ideal_o) / float(n);
        float dash    = arc / float(DASHED_PER_BOND) * RING_DASH_FILL;
        dashes = dash >= RING_MIN_DASH * 2.0 * r ? DASHED_PER_BOND : 1u;
    }

    // Frame of the torus: x towards the first atom, z along the normal
    vec3 q0  = (p[0] - c) - nrm * dot(p[0] - c, nrm);
    vec3 tng = dot(q0, q0) > 0.0 ? normalize(q0) : normalize(cross(nrm, abs(nrm.x) < 0.9 ? vec3(1.0, 0.0, 0.0) : vec3(0.0, 1.0, 0.0)));
    vec3 bit = cross(nrm, tng);

    vec4 angle[2] = vec4[2](vec4(0.0), vec4(0.0));
    for (int k = 1; k < n; ++k) {
        vec3  q = p[k] - c;
        float a = atan(dot(q, bit), dot(q, tng));
        angle[k / 4][k % 4] = a < 0.0 ? a + 2.0 * PI : a;
    }

    out_frag.center  = c;
    out_frag.normal  = nrm;
    out_frag.tangent = tng;
    out_frag.radii   = vec2(R, r);
    out_frag.angle   = angle;
    out_frag.color   = uvec4[2](uvec4(col[0], col[1], col[2], col[3]), uvec4(col[4], col[5], col[6], 0u));
    out_frag.count_picking = uvec2(uint(n), u_ring_base_index + ring_idx);
    out_frag.view_vel = mat3(u_world_to_view) * instance_velocity(local_center / float(n), vel / float(n));
    out_frag.flags    = common_flags;
    out_frag.dashes   = dashes;

    // The box of the torus in view space: half extents R + r along the tangent and bitangent, r along the normal
    vec3 ext = abs(tng) * (R + r) + abs(bit) * (R + r) + abs(nrm) * r;
    vec3 box_lo = c - ext;
    vec3 box_hi = c + ext;

    // Distance from the eye to the near plane (z = -near)
#if ORTHO
    float near = (u_view_to_clip[3][2] + 1.0) / u_view_to_clip[2][2];
#else
    float near = u_view_to_clip[3][2] / (u_view_to_clip[2][2] - 1.0);
#endif

    if (box_lo.z >= -near) {
        // Entirely in front of the near plane
        cull();
        return;
    }
    // The point of the box closest to the image plane is in front of the torus
    float z_front = box_hi.z;
    float z_quad  = min(z_front, -near);

    // The sides of the quad are moved out by a fraction of a pixel, as in bond.vert
    vec2 px = 2.0 / (u_viewport_size * vec2(u_view_to_clip[0][0], u_view_to_clip[1][1]));
    float margin = u_viewport_size.x > 0.0 ? SUBPIXEL_MARGIN * max(abs(px.x), abs(px.y)) : 0.0;

    vec2 img_lo, img_hi;
#if ORTHO
    img_lo = box_lo.xy;
    img_hi = box_hi.xy;
#else
    if (z_front > -near) {
        // Crosses the near plane: the part behind it, within the box in view space
        vec2 sx = box_support(box_lo, box_hi, near, vec2(1.0, 0.0));
        vec2 sy = box_support(box_lo, box_hi, near, vec2(0.0, 1.0));
        img_lo = vec2(sx.x, sy.x);
        img_hi = vec2(sx.y, sy.y);
    } else {
        // The corners of the oriented box, all behind the near plane: their images bound that of the box
        img_lo = vec2( 3.4e38);
        img_hi = vec2(-3.4e38);
        for (int i = 0; i < 8; ++i) {
            vec3 v = c + tng * (R + r) * ((i & 1) != 0 ? 1.0 : -1.0)
                       + bit * (R + r) * ((i & 2) != 0 ? 1.0 : -1.0)
                       + nrm * r       * ((i & 4) != 0 ? 1.0 : -1.0);
            vec2 s = v.xy / -v.z;
            img_lo = min(img_lo, s);
            img_hi = max(img_hi, s);
        }
    }
#endif
    img_lo -= margin;
    img_hi += margin;

    bool hi = (corner & 1u) != 0u;
    bool up = (corner & 2u) != 0u;
    vec2 st = vec2(hi ? img_hi.x : img_lo.x, up ? img_hi.y : img_lo.y);

#if ORTHO
    vec3 view_pos = vec3(st, z_quad);
#else
    vec3 view_pos = vec3(st, -1.0) * -z_quad;
#endif

    vec4 clip = u_view_to_clip * vec4(view_pos, 1.0);
    // Keep the quad behind the near plane, the depth does not change its image
    clip.z = max(clip.z, -clip.w * (1.0 - 1.0e-6));

    out_frag.view_pos = view_pos;
    gl_Position = clip;
}
