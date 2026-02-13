#version 330 core
#extension GL_ARB_shading_language_packing : enable

#ifndef NUM_SUBDIVISIONS
#define NUM_SUBDIVISIONS 8
#endif

#ifndef ORTHONORMALIZE
#define ORTHONORMALIZE 1
#endif

#ifndef USE_B_SPLINE
#define USE_B_SPLINE 1
#endif

uniform float u_tension = 0.5;

uniform usamplerBuffer u_buf_control_point_words;                   // Interleaved gl_control_point_t as uint32 words
uniform usamplerBuffer u_buf_subdivision_cp_indices;                // uvec4(cp0, cp1, cp2, cp3)
uniform samplerBuffer  u_buf_subdivision_params;                    // vec4(t, _, _, _)

out vec3  out_position;
out uint  out_atom_idx;
out vec3  out_velocity;
out float out_segment_t;
out uint  out_secondary_structure_and_flags;
out uvec3 out_support_and_tangent_vector;

vec3 catmull_rom(in vec3 p0, in vec3 p1, in vec3 p2, in vec3 p3, float s, float tension) {
    vec3 v0 = (p2 - p0) * tension;
    vec3 v1 = (p3 - p1) * tension;

    vec3 a = 2.0 * (p1 - p2) + (v0 + v1);
    vec3 b = 3.0 * (p2 - p1) - (2.0 * v0 + v1);

    vec3 res_0 = (a * s * s * s) + (b * s * s);
    vec3 res_1 = (v0 * s) + p1;
    return res_0 + res_1;
}

vec3 catmull_rom_tangent(in vec3 p0, in vec3 p1, in vec3 p2, in vec3 p3, float s, float tension) {
    vec3 v0 = (p2 - p0) * tension;
    vec3 v1 = (p3 - p1) * tension;
    return (2.0 * p1 - 2.0 * p2 + v0 + v1) * 3.0 * s * s + (-3.0 * p1 + 3.0 * p2 - 2.0 * v0 - v1) * 2.0 * s + v0;
}

float b3(float t) {
    float at = abs(t);
    float at1 = 1.0 - at;
    float at2 = 2.0 - at;
    float t1 = at1 * at1 * at1 * (2.0 / 3.0);
    float t2 = at2 * at2 * at2 * (1.0 / 6.0);
    return mix(t2 - t1, t2, step(1.0, at));
}

float b3_t(float t) {
    float at = abs(t);
    float at1 = 1.0 - at;
    float at2 = 2.0 - at;
    float t1 = -sign(t) * at1 * at1 * 2.0;
    float t2 = -sign(t) * at2 * at2 * 0.5;
    return mix(t2 - t1, t2, step(1.0, at));
}

vec3 b_spline(in vec3 p0, in vec3 p1, in vec3 p2, in vec3 p3, float s) {
    return  p0 * b3(s + 1.0) +
            p1 * b3(s) +
            p2 * b3(s - 1.0) +
            p3 * b3(s - 2.0);
}

vec3 b_spline_tangent(in vec3 p0, in vec3 p1, in vec3 p2, in vec3 p3, float s) {
    return  p0 * b3_t(s + 1.0) +
            p1 * b3_t(s) +
            p2 * b3_t(s - 1.0) +
            p3 * b3_t(s - 2.0);
}

vec3 spline(in vec3 p0, in vec3 p1, in vec3 p2, in vec3 p3, float s) {
#if USE_B_SPLINE
    return b_spline(p0, p1, p2, p3, s);
#else
    return catmull_rom(p0, p1, p2, p3, s, u_tension);
#endif
}

vec3 spline_tangent(in vec3 p0, in vec3 p1, in vec3 p2, in vec3 p3, float s) {
#if USE_B_SPLINE
    return b_spline_tangent(p0, p1, p2, p3, s);
#else
    return catmull_rom_tangent(p0, p1, p2, p3, s, u_tension);
#endif
}

vec3 safe_normalize(vec3 v, vec3 fallback) {
    float d2 = dot(v, v);
    return d2 > 0.0 ? v * inversesqrt(d2) : fallback;
}

vec3 unpack_support(uint w0, uint w1, uint w2) {
    vec2 xy = unpackSnorm2x16(w0);
    vec2 zx = unpackSnorm2x16(w1);
    return vec3(xy.x, xy.y, zx.x);
}

vec3 unpack_tangent(uint w0, uint w1, uint w2) {
    vec2 zx = unpackSnorm2x16(w1);
    vec2 yz = unpackSnorm2x16(w2);
    return vec3(zx.y, yz.x, yz.y);
}

void load_cp(uint cp_idx, out vec3 pos, out uint atom_idx, out vec3 vel, out float seg_t, out vec3 ss, out uint flags, out uvec3 svt) {
    int base = int(cp_idx * 12u);

    pos = vec3(
        uintBitsToFloat(texelFetch(u_buf_control_point_words, base + 0).r),
        uintBitsToFloat(texelFetch(u_buf_control_point_words, base + 1).r),
        uintBitsToFloat(texelFetch(u_buf_control_point_words, base + 2).r)
    );

    atom_idx = texelFetch(u_buf_control_point_words, base + 3).r;

    vel = vec3(
        uintBitsToFloat(texelFetch(u_buf_control_point_words, base + 4).r),
        uintBitsToFloat(texelFetch(u_buf_control_point_words, base + 5).r),
        uintBitsToFloat(texelFetch(u_buf_control_point_words, base + 6).r)
    );

    seg_t = uintBitsToFloat(texelFetch(u_buf_control_point_words, base + 7).r);

    uint ssf = texelFetch(u_buf_control_point_words, base + 8).r;
    ss = unpackUnorm4x8(ssf).xyz;
    flags = (ssf >> 24U) & 0xFFU;

    svt = uvec3(
        texelFetch(u_buf_control_point_words, base + 9).r,
        texelFetch(u_buf_control_point_words, base + 10).r,
        texelFetch(u_buf_control_point_words, base + 11).r
    );
}

void main() {
    uvec4 cp_idx = texelFetch(u_buf_subdivision_cp_indices, gl_VertexID);
    float t = clamp(texelFetch(u_buf_subdivision_params, gl_VertexID).x, 0.0, 1.0);

    vec3 cp0, cp1, cp2, cp3;
    vec3 cv0, cv1, cv2, cv3;
    vec3 ss0, ss1, ss2, ss3;
    uint ai_tmp;
    uint ai0, ai1;
    float seg_t0, seg_t1, seg_t2, seg_t3;
    uint fl_tmp;
    uint fl0, fl1;
    uvec3 svt0, svt1, svt2, svt3;

    load_cp(cp_idx.x, cp0, ai_tmp, cv0, seg_t0, ss0, fl_tmp, svt0);
    load_cp(cp_idx.y, cp1, ai0,    cv1, seg_t1, ss1, fl0,    svt1);
    load_cp(cp_idx.z, cp2, ai1,    cv2, seg_t2, ss2, fl1,    svt2);
    load_cp(cp_idx.w, cp3, ai_tmp, cv3, seg_t3, ss3, fl_tmp, svt3);

    vec3 sv0 = unpack_support(svt0.x, svt0.y, svt0.z);
    vec3 sv1 = unpack_support(svt1.x, svt1.y, svt1.z);
    vec3 sv2 = unpack_support(svt2.x, svt2.y, svt2.z);
    vec3 sv3 = unpack_support(svt3.x, svt3.y, svt3.z);

    sv0 *= sign(dot(sv0, sv1));
    sv2 *= sign(dot(sv1, sv2));
    sv3 *= sign(dot(sv2, sv3));

    vec3 s_vec = safe_normalize(spline(sv0, sv1, sv2, sv3, t), vec3(0.0, 0.0, 1.0));
    vec3 s_tan = safe_normalize(spline_tangent(cp0, cp1, cp2, cp3, t), vec3(1.0, 0.0, 0.0));
    vec3 s_sec = catmull_rom(ss0, ss1, ss2, ss3, t, 0.5);
    vec3 s_pos = spline(cp0, cp1, cp2, cp3, t);
    vec3 s_vel = spline(cv0, cv1, cv2, cv3, t);

#if ORTHONORMALIZE
    s_vec = safe_normalize(s_vec - s_tan * dot(s_vec, s_tan), vec3(0.0, 0.0, 1.0));
#endif

    bool first_sample = t <= (0.5 / float(NUM_SUBDIVISIONS));
    bool last_sample  = t >= (1.0 - 0.5 / float(NUM_SUBDIVISIONS));
    uint s_flags = first_sample ? fl0 : (last_sample ? fl1 : 0U);

    out_position = s_pos;
    out_atom_idx = t < 0.5 ? ai0 : ai1;
    out_velocity = s_vel;
    out_segment_t = seg_t1 + t;
    out_secondary_structure_and_flags = (s_flags << 24U) | (packUnorm4x8(vec4(s_sec, 0.0)) & 0x00FFFFFFU);
    out_support_and_tangent_vector[0] = packSnorm2x16(s_vec.xy);
    out_support_and_tangent_vector[1] = packSnorm2x16(vec2(s_vec.z, s_tan.x));
    out_support_and_tangent_vector[2] = packSnorm2x16(s_tan.yz);
}
