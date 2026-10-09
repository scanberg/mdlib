#version 330 core

// Backbone representations (cartoon and ribbons) by vertex pulling.
//
// One instance is drawn per control point (one per residue), covering the part of the spline that belongs to its
// residue: the spline parameter u in [-0.5, 0.5] around the control point. The instance is a tube of S + 1 rings of
// P profile vertices, plus a cap at each end, triangulated by a static index buffer (built in md_gl.c).
//
// The rings are evaluated once per frame by compute_spline_subdivide.vert, R per control point: ring s of control
// point c lies at u = s / R after c, in the frame of c. With stride = R / S, ring j of the instance of c is ring
// (j + S / 2) * stride of the previous control point for j < S / 2 (brought into the frame of c), and ring
// (j - S / 2) * stride of c otherwise. At the ends of a chain the rings beyond the end collapse onto it.
//
// Adjacent instances share the ring at u = +-0.5, which may be flipped between their frames: the profiles are point
// symmetric and profile vertex k + P / 2 is the exact negation of vertex k, so the shared ring has the same set of
// vertices and the tube is closed without cracks.
//
// The surface normal is analytic: it includes the change of the profile size along the spline (transitions between
// secondary structures) and the twist of the frame, so shading matches the actual surface.
//
// The level of detail is chosen per chain (md_gl.c): a coarser mesh has fewer segments, which use every stride:th
// ring, and fewer profile vertices. Chains do not share rings, so chains of different levels do not crack.

#ifndef REP_RIBBONS
#define REP_RIBBONS 0
#endif

#define FLAG_FLIP_NEXT 16u
#define PI 3.14159265358979323846

// Must match compute_spline_subdivide.vert
#define LENGTH_SCALE 16.0
#define TWIST_SCALE 8.0
#define SS_T_SCALE 4.0

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

    vec4  u_scale;  // Cartoon: (coil, helix, sheet, -) scale. Ribbons: (half width, half thickness, -, -)
    uvec4 u_res;    // (S: segments per residue, P: profile vertex count, O: cap outline vertex count, K: ribbon face subdivisions)
    uvec4 u_rings;  // (R: rings per control point in the ring buffer, stride = R / S, -, -)
};

uniform uint u_instance_offset;                 // Control point of instance 0

uniform usamplerBuffer u_buf_rings;             // Ring records (compute_spline_subdivide.vert), 3 x uvec4, R per control point
uniform usamplerBuffer u_buf_control_points;    // Oriented control points, gl_control_point_t as 3 x uvec4
uniform usamplerBuffer u_buf_neighbors;         // uvec4(prev, next, prev2, next2) per control point, clamped to the chain
uniform samplerBuffer  u_atom_color_buffer;
uniform usamplerBuffer u_atom_flags_buffer;

out Fragment {
    smooth vec3 view_coord;
    smooth vec3 view_velocity;
    smooth vec4 color;
    smooth vec3 view_normal;
    flat   uint picking_idx;
} out_frag;

float unpack_unorm8 (uint v) { return float(v & 0xFFu) * (1.0 / 255.0); }
float unpack_snorm8 (uint v) { return clamp(float(int(v << 24u) >> 24) * (1.0 / 127.0), -1.0, 1.0); }
float unpack_unorm16(uint v) { return float(v & 0xFFFFu) * (1.0 / 65535.0); }
float unpack_snorm16(uint v) { return clamp(float(int(v << 16u) >> 16) * (1.0 / 32767.0), -1.0, 1.0); }

struct Ring {
    vec3  position;
    vec3  velocity;
    vec3  x;        // Support (major axis of the profile)
    vec3  z;        // Tangent
    float len_t;    // |d position / du|
    float twist;    // Rotation of the frame around z per unit u
    vec3  ss;
    vec3  ss_t;
};

Ring load_ring(uint idx) {
    int base = int(idx) * 3;
    uvec4 w0 = texelFetch(u_buf_rings, base + 0);
    uvec4 w1 = texelFetch(u_buf_rings, base + 1);
    uvec4 w2 = texelFetch(u_buf_rings, base + 2);

    Ring r;
    r.position = uintBitsToFloat(w0.xyz);
    r.velocity = uintBitsToFloat(uvec3(w0.w, w1.xy));
    r.len_t    = unpack_unorm16(w1.z) * LENGTH_SCALE;
    r.twist    = unpack_snorm16(w1.z >> 16u) * TWIST_SCALE;
    r.ss_t     = vec3(unpack_snorm8(w1.w), unpack_snorm8(w1.w >> 8u), unpack_snorm8(w1.w >> 16u)) * SS_T_SCALE;
    r.x        = vec3(unpack_snorm16(w2.x), unpack_snorm16(w2.x >> 16u), unpack_snorm16(w2.y));
    r.z        = vec3(unpack_snorm16(w2.y >> 16u), unpack_snorm16(w2.z), unpack_snorm16(w2.z >> 16u));
    r.ss       = vec3(unpack_unorm8(w2.w), unpack_unorm8(w2.w >> 8u), unpack_unorm8(w2.w >> 16u));
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

// Profile vertex in the cross section plane of the frame (x, y = cross(z, x)).
// pos:   position
// dir:   direction of the profile curve, counter clockwise around z (any positive scale)
// pos_t: derivative of the position along the spline
struct ProfileVertex {
    vec2 pos;
    vec2 dir;
    vec2 pos_t;
};

#if REP_RIBBONS
// Box of half extents u_scale.xy, the wide faces are subdivided in K segments.
// Vertices: bottom face K + 1 (-x to +x), right face 2, then the point reflection of these (top face, left face).
// The corners are duplicated for flat normals.
ProfileVertex profile_vertex(uint k, vec2 extent, vec2 extent_t) {
    uint K = u_res.w;
    uint half_count = K + 3u;
    float sgn = k < half_count ? 1.0 : -1.0;
    uint kk = k < half_count ? k : k - half_count;

    vec2 p, d;
    if (kk <= K) {
        p = vec2(-1.0 + 2.0 * float(kk) / float(K), -1.0);
        d = vec2(1.0, 0.0);
    } else {
        p = vec2(1.0, kk == K + 1u ? -1.0 : 1.0);
        d = vec2(0.0, 1.0);
    }

    ProfileVertex v;
    v.pos   = sgn * p * extent;
    v.dir   = sgn * d;
    v.pos_t = sgn * p * extent_t;
    return v;
}

// Outline of the cap, counter clockwise: the bottom face, then the top face
vec2 outline_vertex(uint m, vec2 extent) {
    uint K = u_res.w;
    uint half_count = K + 1u;
    float sgn = m < half_count ? 1.0 : -1.0;
    uint mm = m < half_count ? m : m - half_count;
    return sgn * vec2(-1.0 + 2.0 * float(mm) / float(K), -1.0) * extent;
}

#else
// Ellipse of semi axes extent. The angles are distributed between uniform in angle (round profile) and uniform along
// the major axis (flat profile), which keeps the faces of a flat, twisted ribbon from folding over wide triangles
// while the round edges keep enough vertices.
vec2 ellipse_dir(uint k, vec2 extent) {
    uint P = u_res.y;
    uint half_count = P / 2u;
    float sgn = k < half_count ? 1.0 : -1.0;
    float f = float(k < half_count ? k : k - half_count) / float(half_count);
    float flatness = 1.0 - min(extent.x, extent.y) / max(max(extent.x, extent.y), 1.0e-6);
    float theta = mix(f * PI, acos(1.0 - 2.0 * f), flatness);
    return sgn * vec2(cos(theta), sin(theta));
}

ProfileVertex profile_vertex(uint k, vec2 extent, vec2 extent_t) {
    vec2 cs = ellipse_dir(k, extent);
    ProfileVertex v;
    v.pos   = cs * extent;
    v.dir   = vec2(-cs.y * extent.x, cs.x * extent.y);
    v.pos_t = cs * extent_t;
    return v;
}

vec2 outline_vertex(uint m, vec2 extent) {
    return ellipse_dir(m, extent) * extent;
}

#endif

void emit_hidden() {
    out_frag.view_coord    = vec3(0.0);
    out_frag.view_velocity = vec3(0.0);
    out_frag.color         = vec4(0.0);
    out_frag.view_normal   = vec3(0.0, 0.0, 1.0);
    out_frag.picking_idx   = 0u;
    // Every vertex of a hidden primitive gets the same position outside of the clip volume
    gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
}

void main() {
    uint S = u_res.x;
    uint P = u_res.y;
    uint O = u_res.z;
    uint half_S = S / 2u;

    uint R = u_rings.x;
    uint stride = u_rings.y;

    uint c = u_instance_offset + uint(gl_InstanceID);
    uint atom_idx = load_atom_idx(c);

    uvec4 nb = texelFetch(u_buf_neighbors, int(c));
    bool chain_beg = nb.x == c;
    bool chain_end = nb.y == c;

    if ((chain_beg && chain_end) || !atom_visible(atom_idx)) {
        emit_hidden();
        return;
    }

    // Decode the vertex: tube vertices [0, (S + 1) * P), then the caps (center + outline) at the beginning and the end
    uint vid = uint(gl_VertexID);
    uint tube_count = (S + 1u) * P;
    bool is_cap = vid >= tube_count;
    uint j, k = 0u, cap = 0u, cap_vertex = 0u;
    if (!is_cap) {
        j = vid / P;
        k = vid - j * P;
    } else {
        uint cv = vid - tube_count;
        cap = cv / (O + 1u);
        cap_vertex = cv - cap * (O + 1u);
        j = cap == 0u ? 0u : S;
        // A cap closes the tube at the end of a chain or next to a hidden residue
        bool open = cap == 0u ? (chain_beg || !atom_visible(load_atom_idx(nb.x))) : (chain_end || !atom_visible(load_atom_idx(nb.y)));
        if (!open) {
            emit_hidden();
            return;
        }
    }

    // Ring of the instance
    uint ring_idx;
    float sgn = 1.0;
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
    Ring r = load_ring(ring_idx);
    vec3 x = r.x * sgn;
    vec3 z = r.z;
    vec3 y = cross(z, x);

#if REP_RIBBONS
    vec2 extent   = u_scale.xy;
    vec2 extent_t = vec2(0.0);
#else
    const vec2 coil_scale  = vec2(0.2, 0.2);
    const vec2 helix_scale = vec2(1.2, 0.2);
    const vec2 sheet_scale = vec2(1.5, 0.1);
    vec2 cs = coil_scale  * u_scale.x;
    vec2 hs = helix_scale * u_scale.y;
    vec2 ss = sheet_scale * u_scale.z;
    vec2 extent   = max(r.ss.x * cs + r.ss.y * hs + r.ss.z * ss, vec2(1.0e-3));
    vec2 extent_t = r.ss_t.x * cs + r.ss_t.y * hs + r.ss_t.z * ss;
#endif

    vec3 position;
    vec3 normal;
    if (!is_cap) {
        ProfileVertex pv = profile_vertex(k, extent, extent_t);
        // n = d/dprofile x d/du of the surface p(u) + x(u) pos.x + y(u) pos.y, with the frame rotating by twist around z
        vec3 n_local = vec3(r.len_t * pv.dir.y, -r.len_t * pv.dir.x, pv.dir.x * pv.pos_t.y - pv.dir.y * pv.pos_t.x + r.twist * dot(pv.pos, pv.dir));
        position = r.position + x * pv.pos.x + y * pv.pos.y;
        normal   = x * n_local.x + y * n_local.y + z * n_local.z;
    } else {
        vec2 p = cap_vertex == 0u ? vec2(0.0) : outline_vertex(cap_vertex - 1u, extent);
        position = r.position + x * p.x + y * p.y;
        normal   = cap == 0u ? -z : z;
    }

    vec4 view_coord = u_world_to_view * vec4(position, 1.0);
    out_frag.view_coord    = view_coord.xyz;
    out_frag.view_velocity = vec3(u_world_to_view * vec4(r.velocity, 0.0));
    out_frag.color         = texelFetch(u_atom_color_buffer, int(atom_idx));
    out_frag.view_normal   = mat3(u_world_to_view_normal) * normal;
    // The instance is the part of the spline that belongs to one residue: it picks as that backbone segment (control
    // point c is segment c), not as its CA, which other representations may show and pick as an atom of its own
    out_frag.picking_idx   = u_backbone_base_index + c;
    gl_Position = u_view_to_clip * view_coord;
}
