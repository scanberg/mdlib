#version 330 core

// Backbone representations with the stadium profile (cartoon, tubes) ray cast, after Lindow et al. (Interactive
// Visualization of RNA and DNA Structures, 2019): the alternative to the mesh of backbone.vert, which costs little per
// control point and per pixel what the backbone covers.
//
// The spline is cut into straight segments between its rings (compute_spline_subdivide.vert), and each segment is one
// quad: the instance of control point c covers u in [-0.5, 0.5] about it with up to S segments (S quads, u_res.x), as an
// instance of backbone.vert does with its S + 1 rings (see there: the same rings, flipped into the frame of c, collapsed
// at the ends of a chain). How many it uses depends on how much the spline turns there and how much its profile
// changes (segment_step): straight segments of a curve crease where they meet, by the angle the curve turns over one,
// which the ambient occlusion shows as bands across a helix (it turns by 80 to 90 degrees a residue). A helix gets all
// S, a strand or a straight coil as few as 2, each turning by at most u_scale.w radians. The rest of the quads collapse.
//
// At a ring the profile is a stadium, the core segment from P - l X to P + l X (P the position, X the support) widened
// by r, its flat sides facing +-Y = Z x X. The segment between rings 0 and 1 is the ruled surface between the cores
// thickened by the radius, blended from r0 to r1: backbone_rc.frag ray casts it in closed form (see there). A circle
// (l = 0) makes it one round cone, exactly. Where the band ends a segment is closed: round at the end of a chain, flat
// (cut at the plane of the ring) where the neighbor is not drawn (its control point hidden, or not of this pass: a pass
// of the highlighted atoms only must give the same depth where it draws as the pass of all did); elsewhere the next
// segment closes it.
//
// A segment is cut at the planes of its rings (backbone_rc.frag), but where the chain ends, so the segments of a
// chain do not overlap: a quad bounds its own part of the band, and no pixel is ray cast twice where they meet (a
// segment of a helix is 0.3 A long, and its quad reached a radius past either ring; a lipid tail is 2.2 A thick). Between
// the planes it lies within the hull of its cross sections there: the stadium of the ring, whose round ends are the
// round cones of the sides cut by the plane, ellipses a little longer than the circles (bounded here by the rectangle of
// the stadium widened by 1 / cos of the angle between a side and the tangent, and the half angle of the cone where the
// radius changes). Where the band ends, the rung capsule and the end spheres reach out by r past the plane: the box
// out to there. The image of the segment lies within the hull of the images of those corners (a central projection keeps
// convex hulls in front of the eye), and the quad is the rectangle about them along the image of the cores or of the
// axis of the segment, whichever is smaller: a band is a long thin parallelogram on screen at any angle, which a
// rectangle along the axes of the screen bounds poorly. Across the near plane, the box of the corners bounds it instead.
// It lies at the depth of the point of the segment's box closest to the image plane, or at the near plane if that is in
// front of it, as ring.vert places its quads: the ray of a fragment starts on it, in front of the segment.

#ifndef ORTHO
#define ORTHO 0
#endif

#ifndef REP_TUBES
#define REP_TUBES 0
#endif

#define FLAG_FLIP_NEXT 16u

// Must match compute_spline_subdivide.vert
#define SS_T_SCALE 4.0
#define TUBE_RADIUS_MAX 8.0

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

    vec4  u_scale;  // Cartoon: (coil, helix, sheet, max turn) scale. Tubes: (scale, band radius scale, -, max turn). The
                    // largest angle a segment may turn by (radians), 0 for S segments everywhere
    uvec4 u_res;    // (S: most segments per control point, -, -, -)
    uvec4 u_rings;  // (R: rings per control point in the ring buffer, stride = R / S, -, -)
};

uniform uint u_instance_offset;                 // Control point of instance 0
uniform vec2 u_viewport_size;                   // Pixels

uniform usamplerBuffer u_buf_rings;             // Ring records (compute_spline_subdivide.vert), 3 x uvec4, R per control point
uniform usamplerBuffer u_buf_control_points;    // Oriented control points, gl_control_point_t as 3 x uvec4
uniform usamplerBuffer u_buf_neighbors;         // uvec4(prev, next, prev2, next2) per control point, clamped to the chain
uniform samplerBuffer  u_atom_color_buffer;
uniform usamplerBuffer u_atom_flags_buffer;

out Fragment {
    flat vec3  core[4];         // View space: the core of ring 0 (a0, a1), then that of ring 1 (b0, b1)
    flat vec2  radius;          // r0, r1
    flat vec3  normal[2];       // View space normals of the flat sides at the rings (Y), unit, on the same side
    flat vec3  vel[2];          // View space velocities at the rings
    flat vec3  tangent[2];      // View space tangents of the spline at the rings, unit
    flat uvec2 open_end;        // How the band ends at ring 0, at ring 1: 0 not, 1 round (a chain), 2 flat (a neighbor not drawn)
    flat vec4  color;
    flat uint  picking_idx;
    smooth vec3 view_pos;       // On the quad, the ray of the fragment starts here
} out_frag;

float unpack_unorm8 (uint v) { return float(v & 0xFFu) * (1.0 / 255.0); }
float unpack_snorm16(uint v) { return clamp(float(int(v << 16u) >> 16) * (1.0 / 32767.0), -1.0, 1.0); }

struct Ring {
    vec3  position;
    vec3  velocity;
    vec3  x;        // Support
    vec3  z;        // Tangent
    vec3  ss;
    float radius;   // Of a tube
};

Ring load_ring(uint idx) {
    int base = int(idx) * 3;
    uvec4 w0 = texelFetch(u_buf_rings, base + 0);
    uvec4 w1 = texelFetch(u_buf_rings, base + 1);
    uvec4 w2 = texelFetch(u_buf_rings, base + 2);

    Ring r;
    r.position = uintBitsToFloat(w0.xyz);
    r.velocity = uintBitsToFloat(uvec3(w0.w, w1.xy));
    r.x        = vec3(unpack_snorm16(w2.x), unpack_snorm16(w2.x >> 16u), unpack_snorm16(w2.y));
    r.z        = vec3(unpack_snorm16(w2.y >> 16u), unpack_snorm16(w2.z), unpack_snorm16(w2.z >> 16u));
    r.ss       = vec3(unpack_unorm8(w2.w), unpack_unorm8(w2.w >> 8u), unpack_unorm8(w2.w >> 16u));
    r.radius   = unpack_unorm8(w2.w >> 24u) * TUBE_RADIUS_MAX;
    return r;
}

uint load_atom_idx(uint cp_idx) {
    return texelFetch(u_buf_control_points, int(cp_idx) * 3).w;
}

uint load_flags(uint cp_idx) {
    return (texelFetch(u_buf_control_points, int(cp_idx) * 3 + 2).x >> 24u) & 0xFFu;
}

bool atom_visible(uint atom_idx) {
    uint flags = texelFetch(u_atom_flags_buffer, int(atom_idx)).x;
    float alpha = texelFetch(u_atom_color_buffer, int(atom_idx)).a;
    return (flags & u_atom_mask) == u_atom_mask && alpha > 0.0;
}

// (l, r) of the stadium at a ring, as backbone.vert takes them
vec2 stadium(Ring r) {
#if REP_TUBES
    vec2 band = vec2(1.0, r.ss.y > 0.0 ? u_scale.y : 1.0) * u_scale.x;
    return max(vec2(r.ss.y * TUBE_RADIUS_MAX, r.radius) * band, vec2(0.0, 1.0e-3));
#else
    const vec2 coil_scale  = vec2(0.0, 0.2);
    const vec2 helix_scale = vec2(1.0, 0.2);
    const vec2 sheet_scale = vec2(1.4, 0.1);
    return max(r.ss.x * coil_scale * u_scale.x + r.ss.y * helix_scale * u_scale.y + r.ss.z * sheet_scale * u_scale.z, vec2(0.0, 1.0e-3));
#endif
}

// The tangent and the stadium (l, r) of a ring, from the last word of its record
void ring_shape(uint idx, out vec3 z, out vec2 lr) {
    uvec4 w2 = texelFetch(u_buf_rings, int(idx) * 3 + 2);
    Ring r;
    r.position = vec3(0.0);
    r.velocity = vec3(0.0);
    r.x        = vec3(0.0);
    r.z        = vec3(unpack_snorm16(w2.y >> 16u), unpack_snorm16(w2.z), unpack_snorm16(w2.z >> 16u));
    r.ss       = vec3(unpack_unorm8(w2.w), unpack_unorm8(w2.w >> 8u), unpack_unorm8(w2.w >> 16u));
    r.radius   = unpack_unorm8(w2.w >> 24u) * TUBE_RADIUS_MAX;
    z  = normalize(r.z);
    lr = stadium(r);
}

// Ring j in [0, S] of the instance of control point c, and the sign that brings its support into the frame of c
void instance_ring(uint c, uint j, uvec4 nb, bool chain_beg, bool chain_end, out uint ring_idx, out float sgn) {
    uint S = u_res.x;
    uint half_S = S / 2u;
    uint R = u_rings.x;
    uint stride = u_rings.y;
    sgn = 1.0;
    if (j < half_S) {
        if (chain_beg) {
            ring_idx = c * R;
        } else {
            ring_idx = nb.x * R + (j + half_S) * stride;
            sgn = (load_flags(nb.x) & FLAG_FLIP_NEXT) != 0u ? -1.0 : 1.0;
        }
    } else {
        ring_idx = chain_end ? c * R : c * R + (j - half_S) * stride;
    }
}

// The rings a segment of the instance of c spans: the most for which each segment turns by at most u_scale.w and its
// profile changes by at most a tenth of its extent, from the turn and the change over each half of the instance (rings
// 0, S / 2 and S). It divides S / 2, so no segment straddles the control point, nor the ends of a chain.
uint segment_step(uint c, uvec4 nb, bool chain_beg, bool chain_end) {
    uint S = u_res.x;
    uint half_S = S / 2u;
    if (!(u_scale.w > 0.0) || half_S < 2u) return 1u;

    uint i0, ih, i1;
    float sg;
    instance_ring(c, 0u,     nb, chain_beg, chain_end, i0, sg);
    instance_ring(c, half_S, nb, chain_beg, chain_end, ih, sg);
    instance_ring(c, S,      nb, chain_beg, chain_end, i1, sg);
    vec3 z0, zh, z1;
    vec2 s0, sh, s1;
    ring_shape(i0, z0, s0);
    ring_shape(ih, zh, sh);
    ring_shape(i1, z1, s1);

    float turn   = acos(clamp(dot(z0, zh), -1.0, 1.0)) + acos(clamp(dot(zh, z1), -1.0, 1.0));
    float extent = max(max(s0.x + s0.y, sh.x + sh.y), s1.x + s1.y);
    float change = (length(s0 - sh) + length(sh - s1)) / max(extent, 1.0e-6);
    for (uint k = half_S; k > 1u; --k) {
        if (half_S % k != 0u) continue;
        float n = float(S / k);
        if (turn <= u_scale.w * n && change <= 0.1 * n) return k;
    }
    return 1u;
}

// Support of the image of the part of the box [lo, hi] (view space) behind the near plane, along the unit direction n
// of the image plane z = -1, as (s_min, s_max). As in bond.vert.
vec2 box_support(vec3 lo, vec3 hi, float near, vec2 n) {
    float m  = dot(n, 0.5 * (lo.xy + hi.xy));
    float e  = dot(abs(n), 0.5 * (hi.xy - lo.xy));
    float d0 = near;
    float d1 = max(-lo.z, near);
    return vec2(min((m - e) / d0, (m - e) / d1), max((m + e) / d0, (m + e) / d1));
}

// The image of a point in front of the eye, on the image plane z = -1 (or the view plane)
vec2 image(vec3 p) {
#if ORTHO
    return p.xy;
#else
    return p.xy / -p.z;
#endif
}

void cull() {
    // All corners of the quad coincide: no area, no fragments
    gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
}

void main() {
    uint S      = u_res.x;
    uint half_S = S / 2u;
    uint c      = u_instance_offset + uint(gl_InstanceID);
    uint j      = uint(gl_VertexID) >> 2u;
    uint corner = uint(gl_VertexID) & 3u;

    uint atom_idx = load_atom_idx(c);
    uvec4 nb = texelFetch(u_buf_neighbors, int(c));
    bool chain_beg = nb.x == c;
    bool chain_end = nb.y == c;

    // A chain of one control point has no spline; beyond the ends of a chain the rings collapse: no segment there
    if ((chain_beg && chain_end) || (chain_beg && j < half_S) || (chain_end && j >= half_S) || !atom_visible(atom_idx)) {
        cull();
        return;
    }

    // The segments this instance uses: a quad of every step-th ring, the others collapse
    uint step = segment_step(c, nb, chain_beg, chain_end);
    if (j % step != 0u) {
        cull();
        return;
    }

    uint  idx0, idx1;
    float sgn0, sgn1;
    instance_ring(c, j,        nb, chain_beg, chain_end, idx0, sgn0);
    instance_ring(c, j + step, nb, chain_beg, chain_end, idx1, sgn1);
    Ring r0 = load_ring(idx0);
    Ring r1 = load_ring(idx1);

    vec2 s0 = stadium(r0);
    vec2 s1 = stadium(r1);

    mat3 to_view = mat3(u_world_to_view);
    vec3 p0 = vec3(u_world_to_view * vec4(r0.position, 1.0));
    vec3 p1 = vec3(u_world_to_view * vec4(r1.position, 1.0));
    vec3 x0 = to_view * (r0.x * sgn0);
    vec3 x1 = to_view * (r1.x * sgn1);
    // The supports in the same sense at both rings, so the ruled surface between the cores does not twist through itself
    if (dot(x0, x1) < 0.0) x1 = -x1;
    vec3 z0 = normalize(to_view * r0.z);
    vec3 z1 = normalize(to_view * r1.z);
    vec3 e0 = x0 * s0.x;
    vec3 e1 = x1 * s1.x;

    vec3 a0 = p0 - e0, a1 = p0 + e0;
    vec3 b0 = p1 - e1, b1 = p1 + e1;
    vec3 y0 = normalize(cross(z0, x0));
    vec3 y1 = normalize(cross(z1, x1));

    // The band ends round where the chain does, flat where the segment of the neighboring instance is not drawn
    uint end0 = (chain_beg && j == half_S) ? 1u : (j == 0u && !atom_visible(load_atom_idx(nb.x))) ? 2u : 0u;
    uint end1 = (chain_end && j + step == half_S) ? 1u : (j + step == S && !atom_visible(load_atom_idx(nb.y))) ? 2u : 0u;

    // The corners of the rectangles of the cross sections at the rings, and of the boxes out past them where the band
    // ends round. k widens the round ends to the cuts of the round cones of the sides (1 / cos(angle + half angle of the cone)).
    vec3  d0 = b0 - a0;
    vec3  d1 = b1 - a1;
    float len = max(min(length(d0), length(d1)), 1.0e-6);
    vec3  u0 = d0 / max(length(d0), 1.0e-6);
    vec3  u1 = d1 / max(length(d1), 1.0e-6);
    float cos_side = min(min(abs(dot(u0, z0)), abs(dot(u1, z0))), min(abs(dot(u0, z1)), abs(dot(u1, z1))));
    float k = 1.01 / max(cos_side - abs(s0.y - s1.y) / len, 0.5);
    vec3 corner_pt[16];
    int  n_pt = 0;
    for (int i = 0; i < 2; ++i) {
        vec3  P   = i == 0 ? p0 : p1;
        vec3  X   = (i == 0 ? x0 : x1) * ((i == 0 ? s0.x : s1.x) + k * (i == 0 ? s0.y : s1.y));
        vec3  Y   = (i == 0 ? y0 : y1) * (k * (i == 0 ? s0.y : s1.y));
        bool  end = (i == 0 ? end0 : end1) == 1u;
        vec3  out_dir = i == 0 ? -z0 * s0.y : z1 * s1.y;
        for (int q = 0; q < 4; ++q) {
            vec3 v = P + X * ((q & 1) != 0 ? 1.0 : -1.0) + Y * ((q & 2) != 0 ? 1.0 : -1.0);
            corner_pt[n_pt++] = v;
            if (end) corner_pt[n_pt++] = v + out_dir;
        }
    }
    vec3 box_lo = vec3( 3.4e38);
    vec3 box_hi = vec3(-3.4e38);
    for (int i = 0; i < 16; ++i) {
        if (i >= n_pt) break;
        box_lo = min(box_lo, corner_pt[i]);
        box_hi = max(box_hi, corner_pt[i]);
    }

    out_frag.core        = vec3[4](a0, a1, b0, b1);
    out_frag.radius      = vec2(s0.y, s1.y);
    out_frag.normal      = vec3[2](y0, y1);
    out_frag.vel         = vec3[2](to_view * r0.velocity, to_view * r1.velocity);
    out_frag.tangent     = vec3[2](z0, z1);
    out_frag.open_end    = uvec2(end0, end1);
    out_frag.color       = texelFetch(u_atom_color_buffer, int(atom_idx));
    // As the mesh: the part of the spline that belongs to the control point picks as it
    out_frag.picking_idx = u_backbone_base_index + c;

    // Distance from the eye to the near plane (z = -near)
#if ORTHO
    float near = (u_view_to_clip[3][2] + 1.0) / u_view_to_clip[2][2];
#else
    float near = u_view_to_clip[3][2] / (u_view_to_clip[2][2] - 1.0);
#endif

    if (box_lo.z >= -near) {
        cull();
        return;
    }
    float z_front = box_hi.z;
    float z_quad  = min(z_front, -near);

    vec2 px = 2.0 / (u_viewport_size * vec2(u_view_to_clip[0][0], u_view_to_clip[1][1]));
    float margin = u_viewport_size.x > 0.0 ? SUBPIXEL_MARGIN * max(abs(px.x), abs(px.y)) : 0.0;

    // The rectangle (g, gp) of the image: (sg, st) its extents along g and gp
    vec2 g  = vec2(1.0, 0.0);
    vec2 gp = vec2(0.0, 1.0);
    vec2 sg, st;
#if !ORTHO
    if (z_front > -near) {
        // Across the near plane the segment may reach behind the eye, where its image is unbounded: the box of its part
        // behind the near plane, along the axes of the screen
        sg = box_support(box_lo, box_hi, near, g);
        st = box_support(box_lo, box_hi, near, gp);
    } else
#endif
    {
        vec2 im[16];
        for (int i = 0; i < 16; ++i) {
            im[i] = i < n_pt ? image(corner_pt[i]) : vec2(0.0);
        }
        // Along the image of the support (across the band) or of the axis, whichever bounds it the tighter
        vec2 dir[3] = vec2[3](image(p0 + x0) - image(p0) + image(p1 + x1) - image(p1), image(p1) - image(p0), vec2(1.0, 0.0));
        float best = 3.4e38;
        for (int d = 0; d < 3; ++d) {
            float d2 = dot(dir[d], dir[d]);
            if (!(d2 > 1.0e-16)) continue;
            vec2 u = dir[d] * inversesqrt(d2);
            vec2 v = vec2(-u.y, u.x);
            vec2 su = vec2(3.4e38, -3.4e38);
            vec2 sv = vec2(3.4e38, -3.4e38);
            for (int i = 0; i < 16; ++i) {
                if (i >= n_pt) break;
                float a = dot(im[i], u);
                float b = dot(im[i], v);
                su = vec2(min(su.x, a), max(su.y, a));
                sv = vec2(min(sv.x, b), max(sv.y, b));
            }
            float area = (su.y - su.x) * (sv.y - sv.x);
            if (area < best) {
                best = area;
                g = u; gp = v; sg = su; st = sv;
            }
        }
    }
    sg += vec2(-margin, margin);
    st += vec2(-margin, margin);

    bool hi = (corner & 1u) != 0u;
    bool up = (corner & 2u) != 0u;
    vec2 img = (hi ? sg.y : sg.x) * g + (up ? st.y : st.x) * gp;

#if ORTHO
    vec3 view_pos = vec3(img, z_quad);
#else
    vec3 view_pos = vec3(img, -1.0) * -z_quad;
#endif

    vec4 clip = u_view_to_clip * vec4(view_pos, 1.0);
    clip.z = max(clip.z, -clip.w * (1.0 - 1.0e-6));

    out_frag.view_pos = view_pos;
    gl_Position = clip;
}
