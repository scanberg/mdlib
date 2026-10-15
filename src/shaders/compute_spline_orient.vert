#version 330 core
#extension GL_ARB_shading_language_packing : enable

// Orients the backbone control points.
//
// The support vector of control point i is the binormal of the spline's own Frenet frame at that point.
// The subdivision evaluates a uniform cubic B-spline through the CA positions, whose first and second
// derivatives at control point i are (CA[i+1] - CA[i-1]) / 2 and CA[i-1] - 2 CA[i] + CA[i+1]. Hence
//
//     T = CA[i+1] - CA[i-1],   N = CA[i-1] - 2 CA[i] + CA[i+1],   S = normalize(cross(N, T))
//
// S is perpendicular to the spline tangent at the control point by construction, so the
// orthonormalization in the subdivision can never collapse it. The CA virtual bond angle stays well
// below 180 degrees in proteins, so N does not vanish. The carbonyl direction (O - C) from the
// extraction pass is only used as a fallback where N degenerates (very short chains).
//
// S is defined up to sign (it is a director). Whether S[i+1] is used flipped relative to S[i] is
// decided here, once per pair, with temporal hysteresis: the decision from the previous frame is kept
// until the angle under it exceeds 90 degrees by a margin. Re-deciding it from scratch every frame
// (nearest direction) makes pairs close to 90 degrees apart flip back and forth with thermal noise,
// which shows up as segments that suddenly rotate during playback. The decision is stored in
// FLAG_FLIP_NEXT and consumed by the subdivision.
//
// A control point with FLAG_FIXED_FACE (the tubes of md_gl_mol_set_tubes with a particle to face, the backbone of a
// nucleic acid facing its bases) takes no orientation from the curve: the flat faces of its profile face the direction
// extracted for it, so the support (the long axis of the profile) is perpendicular to that and to the tangent.

//
// A control point with FLAG_THROUGH is moved so that the B-spline passes (nearly) through it rather than near it: the
// subdivision evaluates (P[i-1] + 4 P[i] + P[i+1]) / 6 at the control points, which shrinks a curve (a helix by 6% of
// its radius at 36 degrees a step, more than half an Ångström for a nucleic acid). Two Jacobi steps towards the P that
// make this the particle x, from P = x, bring that to about a tenth. They are linear, so the velocity takes them as well.

#define FLAG_FLIP_NEXT  16u
#define FLAG_FIXED_FACE 32u
#define FLAG_THROUGH    64u

uniform usamplerBuffer u_buf_control_point_words;       // Extracted control points (gl_control_point_t as uint32 words)
uniform usamplerBuffer u_buf_prev_control_point_words;  // Oriented control points of the previous computation
uniform usamplerBuffer u_buf_neighbors;                 // uvec4(prev, next, prev2, next2), clamped to the chain
uniform int   u_has_history = 0;                        // u_buf_prev_control_point_words holds valid data for this molecule
uniform float u_flip_threshold = -0.64278760968;        // -sin(hysteresis angle), here 40 degrees

layout (location = 0) in vec3  in_position;
layout (location = 1) in uint  in_atom_idx;
layout (location = 2) in vec3  in_velocity;
layout (location = 3) in float in_radius;
layout (location = 4) in uint  in_secondary_structure_and_flags;
layout (location = 5) in uvec3 in_support_and_tangent_vector;

out vec3  out_position;
out uint  out_atom_idx;
out vec3  out_velocity;
out float out_radius;
out uint  out_secondary_structure_and_flags;
out uvec3 out_support_and_tangent_vector;

#ifndef GL_ARB_shading_language_packing
uint packSnorm2x16(in vec2 v) {
    ivec2 iv = ivec2(round(clamp(v, -1.0f, 1.0f) * 32767.0f));
    return uint(iv.y << 16) | uint(iv.x & 0xFFFF);
}

vec2 unpackSnorm2x16(uint p) {
    ivec2 iv = ivec2(int(p << 16) >> 16, int(p) >> 16);
    return clamp(vec2(iv) * (1.0f / 32767.0f), -1.0f, 1.0f);
}
#endif

vec3 load_position(uint cp_idx) {
    int base = int(cp_idx * 12u);
    return vec3(
        uintBitsToFloat(texelFetch(u_buf_control_point_words, base + 0).r),
        uintBitsToFloat(texelFetch(u_buf_control_point_words, base + 1).r),
        uintBitsToFloat(texelFetch(u_buf_control_point_words, base + 2).r)
    );
}

vec3 load_velocity(uint cp_idx) {
    int base = int(cp_idx * 12u);
    return vec3(
        uintBitsToFloat(texelFetch(u_buf_control_point_words, base + 4).r),
        uintBitsToFloat(texelFetch(u_buf_control_point_words, base + 5).r),
        uintBitsToFloat(texelFetch(u_buf_control_point_words, base + 6).r)
    );
}

// Normalized O - C written by the extraction pass
vec3 load_carbonyl(uint cp_idx) {
    int base = int(cp_idx * 12u);
    vec2 xy = unpackSnorm2x16(texelFetch(u_buf_control_point_words, base +  9).r);
    vec2 zw = unpackSnorm2x16(texelFetch(u_buf_control_point_words, base + 10).r);
    return vec3(xy, zw.x);
}

uint load_flags(uint cp_idx) {
    int base = int(cp_idx * 12u);
    return (texelFetch(u_buf_control_point_words, base + 8).r >> 24u) & 0xFFu;
}

uint load_prev_flags(uint cp_idx) {
    int base = int(cp_idx * 12u);
    return (texelFetch(u_buf_prev_control_point_words, base + 8).r >> 24u) & 0xFFu;
}

vec3 any_perpendicular(vec3 t) {
    vec3 a = abs(t.x) < 0.9 ? vec3(1.0, 0.0, 0.0) : vec3(0.0, 1.0, 0.0);
    return normalize(cross(t, a));
}

void frame_at(uint j, out vec3 tangent, out vec3 support) {
    uvec4 nb = texelFetch(u_buf_neighbors, int(j));

    // Ends of a chain borrow the curvature of their neighbour
    uint a = nb.x, b = j, c = nb.y;
    if (nb.x == j) {
        a = j; b = nb.y; c = nb.w;
    } else if (nb.y == j) {
        a = nb.z; b = nb.x; c = j;
    }

    vec3 T = load_position(nb.y) - load_position(nb.x);
    vec3 N = load_position(a) - 2.0 * load_position(b) + load_position(c);

    float t2 = dot(T, T);
    tangent = t2 > 0.0 ? T * inversesqrt(t2) : vec3(0.0, 0.0, 1.0);

    if ((load_flags(j) & FLAG_FIXED_FACE) != 0u) {
        vec3 s = cross(load_carbonyl(j), tangent);
        float s2 = dot(s, s);
        support = s2 > 1.0e-6 ? s * inversesqrt(s2) : any_perpendicular(tangent);
        return;
    }

    vec3 B = cross(N, tangent);
    float len_b = length(B);    // |component of N perpendicular to T| (Angstrom)

    const float lo = 0.25;
    const float hi = 0.75;
    if (len_b >= hi) {
        support = B / len_b;
        return;
    }

    vec3 oc = load_carbonyl(j);
    oc -= tangent * dot(oc, tangent);
    float oc2 = dot(oc, oc);
    vec3 fallback = oc2 > 1.0e-6 ? oc * inversesqrt(oc2) : any_perpendicular(tangent);
    if (len_b <= lo) {
        support = fallback;
        return;
    }

    vec3 s = B / len_b;
    fallback *= dot(fallback, s) < 0.0 ? -1.0 : 1.0;
    vec3 m = mix(fallback, s, smoothstep(lo, hi, len_b));
    support = normalize(m - tangent * dot(m, tangent));
}

void main() {
    out_position = in_position;
    out_atom_idx = in_atom_idx;
    out_velocity = in_velocity;
    out_radius = in_radius;

    uint cp_idx = uint(gl_VertexID);
    uvec4 nb = texelFetch(u_buf_neighbors, int(cp_idx));

    vec3 tangent, support;
    frame_at(cp_idx, tangent, support);

    uint flags = (in_secondary_structure_and_flags >> 24u) & ~FLAG_FLIP_NEXT;

    if ((flags & FLAG_THROUGH) != 0u) {
        // P1[j] = (6 x[j] - x[j-1] - x[j+1]) / 4 at the neighbours, P2[i] = (6 x[i] - P1[i-1] - P1[i+1]) / 4, with the
        // neighbours clamped to the chain
        vec3 xp  = load_position(nb.x),  xn  = load_position(nb.y);
        vec3 xpp = load_position(nb.z),  xnn = load_position(nb.w);
        vec3 vp  = load_velocity(nb.x),  vn  = load_velocity(nb.y);
        vec3 vpp = load_velocity(nb.z),  vnn = load_velocity(nb.w);
        vec3 p1p = (6.0 * xp - xpp - in_position) * 0.25;
        vec3 p1n = (6.0 * xn - in_position - xnn) * 0.25;
        vec3 v1p = (6.0 * vp - vpp - in_velocity) * 0.25;
        vec3 v1n = (6.0 * vn - in_velocity - vnn) * 0.25;
        out_position = (6.0 * in_position - p1p - p1n) * 0.25;
        out_velocity = (6.0 * in_velocity - v1p - v1n) * 0.25;
    }

    if (nb.y != cp_idx) {
        vec3 next_tangent, next_support;
        frame_at(nb.y, next_tangent, next_support);
        float d = dot(support, next_support);

        bool flip = d < 0.0;
        if (u_has_history != 0) {
            bool prev_flip = (load_prev_flags(cp_idx) & FLAG_FLIP_NEXT) != 0u;
            float sigma = prev_flip ? -1.0 : 1.0;
            flip = (sigma * d < u_flip_threshold) ? !prev_flip : prev_flip;
        }
        if (flip) {
            flags |= FLAG_FLIP_NEXT;
        }
    }

    out_secondary_structure_and_flags = (flags << 24u) | (in_secondary_structure_and_flags & 0x00FFFFFFu);
    out_support_and_tangent_vector[0] = packSnorm2x16(support.xy);
    out_support_and_tangent_vector[1] = packSnorm2x16(vec2(support.z, tangent.x));
    out_support_and_tangent_vector[2] = packSnorm2x16(tangent.yz);
}
