#version 330 core

// Backbone representations (cartoon, ribbons and tubes) by vertex pulling.
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
//
// The profile of the cartoon and of the tubes is a stadium (a 2D capsule): the segment from (-l, 0) to (l, 0) of the
// frame (x = support, y = cross(z, x)) widened by a radius r, two half circles joined by flat sides that face +-y. A
// circle is l = 0. The cartoon blends (l, r) by the secondary structure (coil: a circle, helix and sheet: flat bands of
// the same width and thickness as the ellipses that came before), the tubes have the (l, r) of their control points.
// The ribbons are a box.
//
// Tubes (md_gl_mol_set_tubes) are the backbone of a molecule of particles, each control point a particle: l and r are
// those of the control points (times u_scale.x, and the r of a band, l > 0, times u_scale.y too), blended along the
// spline. Their caps are round: the Minkowski sum of the segment and a half sphere of radius r (a half sphere for a
// circle), D + 1 rings from the rim to the top, which is the segment itself moved out by r.

#ifndef REP_RIBBONS
#define REP_RIBBONS 0
#endif

#ifndef REP_TUBES
#define REP_TUBES 0
#endif

#define FLAG_FLIP_NEXT 16u
#define PI 3.14159265358979323846

// Must match compute_spline_subdivide.vert
#define LENGTH_SCALE 16.0
#define TWIST_SCALE 8.0
#define SS_T_SCALE 4.0
#define TUBE_RADIUS_MAX 8.0
#define RADIUS_T_SCALE 8.0

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

    vec4  u_scale;  // Cartoon: (coil, helix, sheet, -) scale. Ribbons: (half width, half thickness, -, -). Tubes: (radius scale, band radius scale, -, -)
    uvec4 u_res;    // (S: segments per residue, P: profile vertex count, O: cap outline vertex count, K: ribbon face subdivisions or D: dome rings of a tube cap)
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
    float radius;   // Of a tube
    float radius_t;
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
    r.radius   = unpack_unorm8(w2.w >> 24u) * TUBE_RADIUS_MAX;
    r.radius_t = unpack_snorm8(w1.w >> 24u) * RADIUS_T_SCALE;
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
// Stadium of extent = (l, r) (see above), and extent_t = d(l, r)/du. The vertices are spread along the outline by the
// angle it turns through on the half circle, and by a radian per half width (l + r) along the flat side: a circle is
// uniform in angle, a flat band gets most of them on its round edges and a couple across each face, and they move
// continuously as l and r change along the spline. Spread uniformly by length, a sheet (l = 14 r) got one vertex on a
// half circle of P = 16: its edge was a slanted chamfer from one face to the other, shaded as if round, which folded
// along the diagonals of its quads where the strand twists. A face needs few: it is flat across, and its normal changes
// linearly across it. Half k < P / 2 goes from (l, -r) around the half circle about (l, 0) to (l, r) and along the side
// at y = r towards (-l, r), where the point reflection (the other half) begins.
// The derivative of a vertex along the spline is taken at a fixed point of the outline: how the vertices slide along it
// is along dir, which the normal does not see.
ProfileVertex profile_vertex(uint k, vec2 extent, vec2 extent_t) {
    uint  H    = u_res.y / 2u;
    float sgn  = k < H ? 1.0 : -1.0;
    float u    = float(k < H ? k : k - H) / float(H);
    float l    = extent.x;
    float r    = extent.y;
    float side = 2.0 * l / (l + r);     // The flat side, in radians of the parameter
    float s    = u * (PI + side);

    ProfileVertex v;
    if (s < PI) {
        float theta = s - 0.5 * PI;
        vec2  cs    = vec2(cos(theta), sin(theta));
        v.pos   = vec2(l, 0.0) + r * cs;
        v.dir   = vec2(-cs.y, cs.x);
        v.pos_t = vec2(extent_t.x, 0.0) + extent_t.y * cs;
    } else {
        v.pos   = vec2(l - (s - PI) * (l + r), r);
        v.dir   = vec2(-1.0, 0.0);
        v.pos_t = vec2(0.0, extent_t.y);
    }
    v.pos   *= sgn;
    v.dir   *= sgn;
    v.pos_t *= sgn;
    return v;
}

vec2 outline_vertex(uint m, vec2 extent) {
    return profile_vertex(m, extent, vec2(0.0)).pos;
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
#if REP_TUBES
        uint cap_size = (u_res.w + 1u) * O;
#else
        uint cap_size = O + 1u;
#endif
        cap = cv / cap_size;
        cap_vertex = cv - cap * cap_size;
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
#elif REP_TUBES
    // (l, r): the half length of the flat sides in the y of the secondary structure, the radius, that of a band scaled
    vec2 band     = vec2(1.0, r.ss.y > 0.0 ? u_scale.y : 1.0) * u_scale.x;
    vec2 extent   = max(vec2(r.ss.y * TUBE_RADIUS_MAX, r.radius) * band, vec2(0.0, 1.0e-3));
    vec2 extent_t = vec2(r.ss_t.y * TUBE_RADIUS_MAX, r.radius_t) * band;
#else
    // (l, r) of the stadiums of coil, helix and sheet: half widths 0.2, 1.2 and 1.5, half thicknesses 0.2, 0.2 and 0.1
    const vec2 coil_scale  = vec2(0.0, 0.2);
    const vec2 helix_scale = vec2(1.0, 0.2);
    const vec2 sheet_scale = vec2(1.4, 0.1);
    vec2 cs = coil_scale  * u_scale.x;
    vec2 hs = helix_scale * u_scale.y;
    vec2 ss = sheet_scale * u_scale.z;
    vec2 extent   = max(r.ss.x * cs + r.ss.y * hs + r.ss.z * ss, vec2(0.0, 1.0e-3));
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
        vec3 out_dir = cap == 0u ? -z : z;
#if REP_TUBES
        // Ring d of the dome at the angle theta from the rim: each vertex of the rim moves towards the point of the segment
        // of the stadium it is r from, and out by r. The last ring is the segment moved out by r.
        uint  D = u_res.w;
        uint  d = cap_vertex / O;
        ProfileVertex pv = profile_vertex(cap_vertex - d * O, extent, vec2(0.0));
        vec2  n2    = vec2(pv.dir.y, -pv.dir.x);    // Outward, unit
        vec2  core  = pv.pos - extent.y * n2;
        float theta = float(d) / float(D) * (0.5 * PI);
        vec2  p     = core + n2 * (extent.y * cos(theta));
        position = r.position + x * p.x + y * p.y + out_dir * (extent.y * sin(theta));
        normal   = (x * n2.x + y * n2.y) * cos(theta) + out_dir * sin(theta);
#else
        vec2 p = cap_vertex == 0u ? vec2(0.0) : outline_vertex(cap_vertex - 1u, extent);
        position = r.position + x * p.x + y * p.y;
        normal   = out_dir;
#endif
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
