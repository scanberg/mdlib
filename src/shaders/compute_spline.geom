#version 330 core
#extension GL_ARB_shading_language_packing : enable

#ifndef NUM_SUBDIVISIONS
#define NUM_SUBDIVISIONS 8
#endif

#ifndef MAX_VERTICES
#define MAX_VERTICES 9
#endif

#ifndef ORTHONORMALIZE
#define ORTHONORMALIZE 1
#endif

#ifndef USE_B_SPLINE
#define USE_B_SPLINE 1
#endif

uniform float u_tension = 0.5; // @NOTE: Only for catmull rom

layout(lines_adjacency) in;
layout(points, max_vertices = MAX_VERTICES) out;

out vec3  out_position;
out uint  out_atom_idx;
out vec3  out_velocity;
out float out_segment_t;
out uint  out_secondary_structure_and_flags;
out uvec3 out_support_and_tangent_vector; // int16[3] * 2 packed into uint32[3]

in Vertex {
    vec3  position;
    uint  atom_index;
    vec3  velocity;
    float segment_t;
    vec3  secondary_structure;
    uint  flags;
    vec3  support_vector;
} in_vert[];

#ifndef GL_ARB_shading_language_packing
uint packSnorm2x16(in vec2 v) {
    ivec2 iv = ivec2(round(clamp(v, -1.0f, 1.0f) * 32767.0f));
    return uint(iv.y << 16) | uint(iv.x & 0xFFFF);
}

uint packUnorm4x8(in vec4 v) {
    uvec4 iv = uvec4(round(clamp(v, 0.0f, 1.0f) * 255.0f));
    return (iv.w << 24U) | (iv.z << 16U) | (iv.y << 8U) | (iv.x);
}
#endif

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

// Inspired by this:
// https://www.it.uu.se/edu/course/homepage/grafik1/ht07/examples/curves.cpp
float b3(float t) {
    float at = abs(t);
    float at1 = 1.0 - at;
    float at2 = 2.0 - at;
    float t1 = at1 * at1 * at1 * (2.0 / 3.0);
    float t2 = at2 * at2 * at2 * (1.0 / 6.0);
    return mix(t2 - t1, t2, step(1.0, at));
}

// b-spline tangent -> 1st derivative
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

vec3 slerp(vec3 v0, vec3 v1, float t) {
    float a = acos(dot(v0/length(v0), v1/length(v1)));
    float iso = 1.0 / sin(a);
    return sin((1.0-t)*a) * iso * v0 + sin(t*a) * iso * v1;
}

void main() {
    vec3 cp[4];
    vec3 cv[4];
    vec3 sv[4];
    vec3 ss[4];
    uint ai[2];
    uint fl[2];

    cp[0] = in_vert[0].position;
    cp[1] = in_vert[1].position;
    cp[2] = in_vert[2].position;
    cp[3] = in_vert[3].position;

    cv[0] = in_vert[0].velocity;
    cv[1] = in_vert[1].velocity;
    cv[2] = in_vert[2].velocity;
    cv[3] = in_vert[3].velocity;

    sv[0] = in_vert[0].support_vector;
    sv[1] = in_vert[1].support_vector;
    sv[2] = in_vert[2].support_vector;
    sv[3] = in_vert[3].support_vector;

    // @NOTE: Flip signs of support vectors with respect to sv[1], this is our reference
    sv[0] *= sign(dot(sv[0], sv[1]));
    sv[2] *= sign(dot(sv[1], sv[2]));
    sv[3] *= sign(dot(sv[2], sv[3]));

    ss[0] = vec3(in_vert[0].secondary_structure);
    ss[1] = vec3(in_vert[1].secondary_structure);
    ss[2] = vec3(in_vert[2].secondary_structure);
    ss[3] = vec3(in_vert[3].secondary_structure);

    ai[0] = in_vert[1].atom_index;
    ai[1] = in_vert[2].atom_index;

    fl[0] = in_vert[1].flags;
    fl[1] = in_vert[2].flags;

    int end_idx = (fl[1] & 2u) != 0U ? NUM_SUBDIVISIONS + 1 : NUM_SUBDIVISIONS;
    int last_idx = end_idx - 1;
    for (int i = 0; i < end_idx; i++) {
        float t = float(i) / float(NUM_SUBDIVISIONS);
        vec3 s_vec = normalize(spline(sv[0], sv[1], sv[2], sv[3], t));
        vec3 s_tan = normalize(spline_tangent(cp[0], cp[1], cp[2], cp[3], t));
        vec3 s_sec = catmull_rom(ss[0], ss[1], ss[2], ss[3], t, 0.5);
        vec3 s_pos = spline(cp[0], cp[1], cp[2], cp[3], t);
        vec3 s_vel = spline(cv[0], cv[1], cv[2], cv[3], t);
        uint s_flags = (i == 0) ? fl[0] : (i == last_idx) ? fl[1] : 0U;

#if ORTHONORMALIZE
        s_vec = normalize(s_vec - s_tan*dot(s_vec, s_tan));
#endif

        out_position = s_pos;
        out_atom_idx = t < 0.5 ? ai[0] : ai[1];   // Pick closest control point for index
        out_velocity = s_vel;
        out_segment_t = in_vert[1].segment_t + t;
        out_secondary_structure_and_flags = (s_flags << 24U) | (packUnorm4x8(vec4(s_sec, 0)) & 0x00FFFFFFU);
        out_support_and_tangent_vector[0] = packSnorm2x16(s_vec.xy);
        out_support_and_tangent_vector[1] = packSnorm2x16(vec2(s_vec.z, s_tan.x));
        out_support_and_tangent_vector[2] = packSnorm2x16(s_tan.yz);
        
        EmitVertex();
        EndPrimitive();
    }
}