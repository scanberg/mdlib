#version 330 core

// Evaluates the rings (cross section frames) of the backbone spline, captured with transform feedback and drawn from by
// backbone.vert, so the spline is evaluated once per ring instead of once per vertex.
//
// Drawn as u_segments points per control point (instanced): ring s of control point c lies on the B-spline segment
// (c, next) at t = s / u_segments, its frame aligned to the support vector of c. The last control point of a chain has
// no segment of its own, its ring 0 is the end of the segment (prev, c) and its remaining rings are unused.
//
// The support vectors are aligned using the relations decided by the orient pass (see compute_spline_orient.vert).
//
// Ring record, 3 x uvec4 (u is the spline parameter, one unit per residue):
// [0] position.xyz, velocity.x (f32)
// [1] velocity.yz (f32), |dp/du| (unorm16, x LENGTH_SCALE) | twist (snorm16, x TWIST_SCALE), d(ss)/du (snorm8 x 3, x SS_T_SCALE)
// [2] x.xy, x.z | z.x, z.yz (snorm16, x = support, z = tangent), ss (unorm8 x 3)
//
// twist is the rotation of the frame around the tangent per unit u. ss holds the secondary structure fractions
// (coil, helix, sheet), blended with a Catmull-Rom spline: they control the size of the cartoon profile.

#define FLAG_FLIP_NEXT 16u

// Must match backbone.vert
#define LENGTH_SCALE 16.0
#define TWIST_SCALE 8.0
#define SS_T_SCALE 4.0

uniform usamplerBuffer u_buf_control_points;    // Oriented control points, gl_control_point_t as 3 x uvec4
uniform usamplerBuffer u_buf_neighbors;         // uvec4(prev, next, prev2, next2) per control point, clamped to the chain
uniform int u_segments;

out uvec4 out_ring_0;
out uvec4 out_ring_1;
out uvec4 out_ring_2;

uint pack_unorm8 (float v) { return uint(round(clamp(v, 0.0, 1.0) * 255.0)); }
uint pack_snorm8 (float v) { return uint(int(round(clamp(v, -1.0, 1.0) * 127.0))) & 0xFFu; }
uint pack_unorm16(float v) { return uint(round(clamp(v, 0.0, 1.0) * 65535.0)); }
uint pack_snorm16(float v) { return uint(int(round(clamp(v, -1.0, 1.0) * 32767.0))) & 0xFFFFu; }

uint pack_snorm16x2(float a, float b) { return pack_snorm16(a) | (pack_snorm16(b) << 16u); }

float unpack_unorm8 (uint v) { return float(v & 0xFFu) * (1.0 / 255.0); }
float unpack_snorm16(uint v) { return clamp(float(int(v << 16u) >> 16) * (1.0 / 32767.0), -1.0, 1.0); }

struct ControlPoint {
    vec3 position;
    vec3 velocity;
    vec3 ss;        // Secondary structure fractions (coil, helix, sheet)
    vec3 support;
    uint flags;
};

ControlPoint load_cp(uint idx) {
    int base = int(idx) * 3;
    uvec4 w0 = texelFetch(u_buf_control_points, base + 0);
    uvec4 w1 = texelFetch(u_buf_control_points, base + 1);
    uvec4 w2 = texelFetch(u_buf_control_points, base + 2);

    ControlPoint cp;
    cp.position = uintBitsToFloat(w0.xyz);
    cp.velocity = uintBitsToFloat(w1.xyz);
    cp.ss       = vec3(unpack_unorm8(w2.x), unpack_unorm8(w2.x >> 8u), unpack_unorm8(w2.x >> 16u));
    cp.flags    = (w2.x >> 24u) & 0xFFu;
    cp.support  = vec3(unpack_snorm16(w2.y), unpack_snorm16(w2.y >> 16u), unpack_snorm16(w2.z));
    return cp;
}

float relation(uint flags) {
    return (flags & FLAG_FLIP_NEXT) != 0u ? -1.0 : 1.0;
}

// Uniform cubic B-spline basis on a segment, t in [0, 1]
vec4 b_spline_weights(float t) {
    float t2 = t * t;
    float t3 = t2 * t;
    float s = 1.0 - t;
    return vec4(s * s * s, 3.0 * t3 - 6.0 * t2 + 4.0, -3.0 * t3 + 3.0 * t2 + 3.0 * t + 1.0, t3) * (1.0 / 6.0);
}

vec4 b_spline_derivative_weights(float t) {
    float t2 = t * t;
    float s = 1.0 - t;
    return vec4(-s * s, 3.0 * t2 - 4.0 * t, -3.0 * t2 + 2.0 * t + 1.0, t2) * 0.5;
}

// Catmull-Rom (tension 0.5) basis and its derivative
vec4 catmull_rom_weights(float t) {
    float t2 = t * t;
    float t3 = t2 * t;
    return vec4(-t3 + 2.0 * t2 - t, 3.0 * t3 - 5.0 * t2 + 2.0, -3.0 * t3 + 4.0 * t2 + t, t3 - t2) * 0.5;
}

vec4 catmull_rom_derivative_weights(float t) {
    float t2 = t * t;
    return vec4(-3.0 * t2 + 4.0 * t - 1.0, 9.0 * t2 - 10.0 * t, -9.0 * t2 + 8.0 * t + 1.0, 3.0 * t2 - 2.0 * t) * 0.5;
}

vec3 blend(vec4 w, vec3 p0, vec3 p1, vec3 p2, vec3 p3) {
    return p0 * w.x + p1 * w.y + p2 * w.z + p3 * w.w;
}

vec3 any_perpendicular(vec3 v) {
    vec3 a = abs(v.x) < 0.9 ? vec3(1.0, 0.0, 0.0) : vec3(0.0, 1.0, 0.0);
    return normalize(cross(v, a));
}

void main() {
    uint c = uint(gl_InstanceID);
    uvec4 nb = texelFetch(u_buf_neighbors, int(c));

    // Segment (cp.y, cp.z) with its clamped neighbours cp.x and cp.w
    uvec4 cp;
    float t;
    bool chain_end = nb.y == c;
    if (chain_end) {
        cp = uvec4(nb.z, nb.x, c, c);
        t = 1.0;
    } else {
        cp = uvec4(nb.x, c, nb.y, nb.w);
        t = float(gl_VertexID) / float(u_segments);
    }

    ControlPoint c0 = load_cp(cp.x);
    ControlPoint c1 = load_cp(cp.y);
    ControlPoint c2 = load_cp(cp.z);
    ControlPoint c3 = load_cp(cp.w);

    // Align the support vectors to c1. Clamped indices at the ends of a chain duplicate the end point.
    vec3 s1 = c1.support;
    vec3 s2 = c2.support * relation(c1.flags);
    vec3 s0 = (cp.x == cp.y) ? s1 : c0.support * relation(c0.flags);
    vec3 s3 = (cp.w == cp.z) ? s2 : c3.support * (relation(c1.flags) * relation(c2.flags));

    vec4 w   = b_spline_weights(t);
    vec4 dw  = b_spline_derivative_weights(t);
    vec4 cw  = catmull_rom_weights(t);
    vec4 dcw = catmull_rom_derivative_weights(t);

    vec3 position   = blend(w,  c0.position, c1.position, c2.position, c3.position);
    vec3 position_t = blend(dw, c0.position, c1.position, c2.position, c3.position);
    vec3 support    = blend(w,  s0, s1, s2, s3);
    vec3 support_t  = blend(dw, s0, s1, s2, s3);
    vec3 velocity   = blend(w,  c0.velocity, c1.velocity, c2.velocity, c3.velocity);
    vec3 ss         = blend(cw,  c0.ss, c1.ss, c2.ss, c3.ss);
    vec3 ss_t       = blend(dcw, c0.ss, c1.ss, c2.ss, c3.ss);

    // The end of a chain is evaluated on the segment of the previous control point: bring it into the frame of c
    if (chain_end) {
        float rel = relation(c1.flags);
        support   *= rel;
        support_t *= rel;
    }

    // Frame: z = tangent, x = support orthogonalized against z
    float len_t = length(position_t);
    vec3 z = len_t > 0.0 ? position_t / len_t : vec3(0.0, 0.0, 1.0);
    vec3 sp = support - z * dot(support, z);
    float len_s = length(sp);
    vec3 x = len_s > 1.0e-6 ? sp / len_s : any_perpendicular(z);
    vec3 y = cross(z, x);
    float twist = len_s > 1.0e-6 ? dot(support_t, y) / len_s : 0.0;

    out_ring_0 = uvec4(floatBitsToUint(position), floatBitsToUint(velocity.x));
    out_ring_1 = uvec4(floatBitsToUint(velocity.yz),
                       pack_unorm16(len_t / LENGTH_SCALE) | (pack_snorm16(twist / TWIST_SCALE) << 16u),
                       pack_snorm8(ss_t.x / SS_T_SCALE) | (pack_snorm8(ss_t.y / SS_T_SCALE) << 8u) | (pack_snorm8(ss_t.z / SS_T_SCALE) << 16u));
    out_ring_2 = uvec4(pack_snorm16x2(x.x, x.y),
                       pack_snorm16x2(x.z, z.x),
                       pack_snorm16x2(z.y, z.z),
                       pack_unorm8(ss.x) | (pack_unorm8(ss.y) << 8u) | (pack_unorm8(ss.z) << 16u));
}
