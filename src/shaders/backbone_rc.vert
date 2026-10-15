#version 330 core

// Backbone representations with the stadium profile (cartoon, tubes) ray cast, after Lindow et al. (Interactive
// Visualization of RNA and DNA Structures, 2019): the alternative to the mesh of backbone.vert, which costs little per
// control point and per pixel what the backbone covers.
//
// The spline is cut into straight segments between its rings (compute_spline_subdivide.vert), S per control point, and
// each segment is one quad: the instance of control point c covers u in [-0.5, 0.5] about it with S segments, as an
// instance of backbone.vert does with its S + 1 rings (see there: the same rings, flipped into the frame of c, collapsed
// at the ends of a chain). At a ring the profile is a stadium, the core segment from P - l X to P + l X (P the position,
// X the support) widened by r. The segment between rings 0 and 1 is the union of the capsules between a point of the
// core of ring 0 and the corresponding point of the core of ring 1, of radius blended from r0 to r1 (a ruled surface
// thickened by the radius): backbone_rc.frag sphere traces it. A circle (l = 0) makes it one capsule, exactly. The
// segments of a chain overlap at their rings, round, so the chain is closed and its ends round without caps.
//
// The quad bounds the image of the box of the segment (view space: the box of the four ends of the cores, out by the
// larger radius), at the depth of the point of that box closest to the image plane, or at the near plane if that is in
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

    vec4  u_scale;  // Cartoon: (coil, helix, sheet, -) scale. Tubes: (scale, -, -, -)
    uvec4 u_res;    // (S: segments per control point, -, -, -)
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
    flat vec3  vel[2];          // View space velocities at the rings
    flat vec3  tangent[2];      // View space tangents of the spline at the rings, unit
    flat uvec2 open_end;        // Ring 0 begins the chain, ring 1 ends it
    flat vec3  box_lo;          // View space box of the segment
    flat vec3  box_hi;
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
    return max(vec2(r.ss.y * TUBE_RADIUS_MAX, r.radius) * u_scale.x, vec2(0.0, 1.0e-3));
#else
    const vec2 coil_scale  = vec2(0.0, 0.2);
    const vec2 helix_scale = vec2(1.0, 0.2);
    const vec2 sheet_scale = vec2(1.4, 0.1);
    return max(r.ss.x * coil_scale * u_scale.x + r.ss.y * helix_scale * u_scale.y + r.ss.z * sheet_scale * u_scale.z, vec2(0.0, 1.0e-3));
#endif
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

    uint  idx0, idx1;
    float sgn0, sgn1;
    instance_ring(c, j,      nb, chain_beg, chain_end, idx0, sgn0);
    instance_ring(c, j + 1u, nb, chain_beg, chain_end, idx1, sgn1);
    Ring r0 = load_ring(idx0);
    Ring r1 = load_ring(idx1);

    vec2 s0 = stadium(r0);
    vec2 s1 = stadium(r1);

    mat3 to_view = mat3(u_world_to_view);
    vec3 p0 = vec3(u_world_to_view * vec4(r0.position, 1.0));
    vec3 p1 = vec3(u_world_to_view * vec4(r1.position, 1.0));
    vec3 e0 = to_view * (r0.x * sgn0) * s0.x;
    vec3 e1 = to_view * (r1.x * sgn1) * s1.x;
    // The cores in the same order at both rings, so the ruled surface between them does not twist through itself
    if (dot(e0, e1) < 0.0) e1 = -e1;

    vec3 a0 = p0 - e0, a1 = p0 + e0;
    vec3 b0 = p1 - e1, b1 = p1 + e1;
    // Out by a little more than the radius: a face which touches the box (one parallel to an axis of the view) is then
    // inside it, and the march starts outside the segment (backbone_rc.frag)
    float r_max = max(s0.y, s1.y) * 1.01;
    vec3 box_lo = min(min(a0, a1), min(b0, b1)) - r_max;
    vec3 box_hi = max(max(a0, a1), max(b0, b1)) + r_max;

    out_frag.core        = vec3[4](a0, a1, b0, b1);
    out_frag.radius      = vec2(s0.y, s1.y);
    out_frag.vel         = vec3[2](to_view * r0.velocity, to_view * r1.velocity);
    out_frag.tangent     = vec3[2](normalize(to_view * r0.z), normalize(to_view * r1.z));
    // The first segment of a chain starts at its beginning, the last ends at its end
    out_frag.open_end    = uvec2(chain_beg && j == half_S ? 1u : 0u, chain_end && j + 1u == half_S ? 1u : 0u);
    out_frag.box_lo      = box_lo;
    out_frag.box_hi      = box_hi;
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

    vec2 img_lo, img_hi;
#if ORTHO
    img_lo = box_lo.xy;
    img_hi = box_hi.xy;
#else
    if (z_front > -near) {
        vec2 sx = box_support(box_lo, box_hi, near, vec2(1.0, 0.0));
        vec2 sy = box_support(box_lo, box_hi, near, vec2(0.0, 1.0));
        img_lo = vec2(sx.x, sy.x);
        img_hi = vec2(sx.y, sy.y);
    } else {
        img_lo = vec2( 3.4e38);
        img_hi = vec2(-3.4e38);
        for (int i = 0; i < 8; ++i) {
            vec3 v = vec3((i & 1) != 0 ? box_hi.x : box_lo.x, (i & 2) != 0 ? box_hi.y : box_lo.y, (i & 4) != 0 ? box_hi.z : box_lo.z);
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
    clip.z = max(clip.z, -clip.w * (1.0 - 1.0e-6));

    out_frag.view_pos = view_pos;
    gl_Position = clip;
}
