#version 330 core

// Slabs (md_gl_mol_set_slabs): flat rounded polygons through particles, the rings of a cartoon (the sugars and bases of
// nucleic acids), raycast within quads by slab.frag. Vertex pulling as ring.vert, from the quad index buffer of the
// bonds: primitive p = u_slab_offset + gl_InstanceID * u_quads_per_instance + (gl_VertexID >> 2) is one quad, a slab.
//
// A slab follows its corners in every frame. In view space:
//   - center: the centroid of the corners
//   - normal: Newell's, the sum of the cross products of consecutive corners about the center, which the corners go
//     counter clockwise about. Robust for a puckered ring (a sugar). Two corners have no plane of their own, and any
//     through them will do: the slab is the capsule between them (a stick), the same in every one.
//   - polygon: the corners projected on the plane through the center, in the 2D frame (tangent towards the first corner,
//     bitangent = normal x tangent)
// The solid is the Minkowski sum of that polygon with a sphere of radius u_radius, with its faces moved in to the half
// thickness u_plate: a rim of capsules along the edges and a plate between them. With u_plate = u_radius the faces are
// flat, with less the plate is recessed within the rim.
//
// A slab is culled when one of its corners is hidden (alpha 0, u_atom_mask), when an edge is longer than a bond can be
// (across the periodic boundary), or when it has no plane. Its color is that of its first corner.
//
// The quad bounds the image of the box of the slab and lies at the depth of the point of that box closest to the image
// plane (or at the near plane, if that is in front of it), so it is in front of every visible point of the slab: the
// raycast starts at the quad, and its depth is a lower bound of the depth slab.frag writes (layout (depth_greater)).

#ifndef ORTHO
#define ORTHO 0
#endif

// Pixels the sides of the quad are moved out by
#define SUBPIXEL_MARGIN (1.0 / 32.0)

// Entries of the slabs (md_gl.c, GL_BUFFER_SLAB): two RGBA32UI texels, the corner count, then up to MAX_CORNERS corners
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

    float u_radius;             // Of the rim
    float u_plate;              // Half thickness of the plate, at most u_radius
    float u_max_d2;             // Longest edge, squared
    uint  u_slab_base_index;    // Picking: a slab writes u_slab_base_index + its index
};

uniform uint u_slab_offset;
uniform uint u_quads_per_instance;
uniform vec2 u_viewport_size;           // Pixels

uniform samplerBuffer  u_buf_atom_pos;      // RGB32F
uniform samplerBuffer  u_buf_atom_vel;      // RGB32F
uniform usamplerBuffer u_buf_atom_flags;    // R8UI
uniform samplerBuffer  u_buf_atom_color;    // RGBA8
uniform usamplerBuffer u_buf_slab;          // RGBA32UI, two per slab

out Fragment {
    flat vec3  center;          // View space
    flat vec3  normal;          // Unit, the corners go counter clockwise about it
    flat vec3  tangent;         // Unit, in the plane towards the first corner
    flat vec4  corner[4];       // The corners in the plane, (x, y) of corner 2i in .xy and of 2i + 1 in .zw
    flat uint  count;
    flat vec4  color;
    flat uint  picking;
    flat vec3  view_vel;        // Mean velocity of the corners
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

void main() {
    uint prim   = u_slab_offset + uint(gl_InstanceID) * u_quads_per_instance + (uint(gl_VertexID) >> 2u);
    uint corner = uint(gl_VertexID) & 3u;

    uvec4 e0 = texelFetch(u_buf_slab, int(2u * prim + 0u));
    uvec4 e1 = texelFetch(u_buf_slab, int(2u * prim + 1u));
    int  n   = int(e0.x);
    uint idx[MAX_CORNERS] = uint[MAX_CORNERS](e0.y, e0.z, e0.w, e1.x, e1.y, e1.z, e1.w);

    if (n < 2 || n > MAX_CORNERS) {
        cull();
        return;
    }

    vec3 p[MAX_CORNERS];
    vec3 c   = vec3(0.0);
    vec3 vel = vec3(0.0);
    for (int k = 0; k < MAX_CORNERS; ++k) {
        if (k >= n) continue;
        uint flags = texelFetch(u_buf_atom_flags, int(idx[k])).x;
        vec4 color = texelFetch(u_buf_atom_color, int(idx[k]));
        if ((flags & u_atom_mask) != u_atom_mask || color.a == 0.0) {
            cull();
            return;
        }
        p[k] = vec3(u_world_to_view * vec4(texelFetch(u_buf_atom_pos, int(idx[k])).xyz, 1.0));
        vel += texelFetch(u_buf_atom_vel, int(idx[k])).xyz;
        c   += p[k];
    }
    c /= float(n);

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
    float size2 = dot(p[0] - c, p[0] - c);
    if (n == 2) {
        // A capsule: the plane through its axis and the coordinate axis furthest from it
        vec3 e = p[1] - p[0];
        vec3 a = abs(e);
        N = cross(e, (a.x <= a.y && a.x <= a.z) ? vec3(1.0, 0.0, 0.0) : (a.y <= a.z ? vec3(0.0, 1.0, 0.0) : vec3(0.0, 0.0, 1.0)));
    }
    float N2 = dot(N, N);
    if (!(N2 > 1.0e-8 * size2 * size2)) {
        cull();
        return;
    }
    vec3 nrm = N * inversesqrt(N2);

    vec3 q0  = (p[0] - c) - nrm * dot(p[0] - c, nrm);
    vec3 tng = dot(q0, q0) > 0.0 ? normalize(q0) : normalize(cross(nrm, abs(nrm.x) < 0.9 ? vec3(1.0, 0.0, 0.0) : vec3(0.0, 1.0, 0.0)));
    vec3 bit = cross(nrm, tng);

    // The polygon in the plane, and its extent
    vec2 q[MAX_CORNERS];
    vec2 q_lo = vec2( 3.4e38);
    vec2 q_hi = vec2(-3.4e38);
    for (int k = 0; k < MAX_CORNERS; ++k) {
        q[k] = vec2(0.0);
        if (k >= n) continue;
        vec3 d = p[k] - c;
        q[k] = vec2(dot(d, tng), dot(d, bit));
        q_lo = min(q_lo, q[k]);
        q_hi = max(q_hi, q[k]);
    }

    float r = u_radius;
    out_frag.center   = c;
    out_frag.normal   = nrm;
    out_frag.tangent  = tng;
    out_frag.corner   = vec4[4](vec4(q[0], q[1]), vec4(q[2], q[3]), vec4(q[4], q[5]), vec4(q[6], vec2(0.0)));
    out_frag.count    = uint(n);
    out_frag.color    = texelFetch(u_buf_atom_color, int(idx[0]));
    out_frag.picking  = u_slab_base_index + prim;
    out_frag.view_vel = mat3(u_world_to_view) * (vel / float(n));

    // The box of the slab: the extent of the polygon and r about it in the plane, r along the normal
    q_lo -= r;
    q_hi += r;
    vec3 box_lo = vec3( 3.4e38);
    vec3 box_hi = vec3(-3.4e38);
    vec3 box_corner[8];
    for (int i = 0; i < 8; ++i) {
        vec3 v = c + tng * ((i & 1) != 0 ? q_hi.x : q_lo.x)
                   + bit * ((i & 2) != 0 ? q_hi.y : q_lo.y)
                   + nrm * ((i & 4) != 0 ? r : -r);
        box_corner[i] = v;
        box_lo = min(box_lo, v);
        box_hi = max(box_hi, v);
    }

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
    // The point of the box closest to the image plane is in front of the slab
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
            vec2 s = box_corner[i].xy / -box_corner[i].z;
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
