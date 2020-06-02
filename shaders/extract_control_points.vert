#version 330 core
#extension GL_ARB_shading_language_packing : enable

uniform samplerBuffer u_buf_atom_pos;
uniform samplerBuffer u_buf_atom_prev_pos; 
uniform samplerBuffer u_buf_secondary_structure;

layout (location = 0) in uint in_residue_idx;
layout (location = 1) in uint in_segment_idx;
layout (location = 2) in uint in_atom_ca;
layout (location = 3) in uint in_atom_c;
layout (location = 4) in uint in_atom_o;

out vec3  out_position;
out uint  out_atom_idx;
out vec3  out_velocity;
out float out_segment_t;
out uint  out_secondary_structure_and_flags;
out uvec3 out_support_and_tangent_vector; // int16[3] * 2 packed into uint32[3]

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

void main() {
    vec3 ca   = texelFetch(u_buf_atom_pos, 		int(in_atom_ca)).xyz;
    vec3 ca_p = texelFetch(u_buf_atom_prev_pos, int(in_atom_ca)).xyz;
    vec3 c    = texelFetch(u_buf_atom_pos, 		int(in_atom_c)).xyz;
    vec3 o    = texelFetch(u_buf_atom_pos, 		int(in_atom_o)).xyz;
    vec4 ss   = texelFetch(u_buf_secondary_structure, int(in_residue_idx));

    vec3 v = normalize(o - c);
    vec3 t = vec3(0);   // Placeholder, tangent is computed analytically in spline shader
    uint f = 0U; // flags

    out_position = ca;
    out_atom_idx = in_atom_ca;
    out_velocity = ca - ca_p;
    out_segment_t = float(in_segment_idx);
    out_secondary_structure_and_flags = (f << 24U) | (packUnorm4x8(ss) & 0x00FFFFFFU);
    out_support_and_tangent_vector[0] = packSnorm2x16(v.xy);
    out_support_and_tangent_vector[1] = packSnorm2x16(vec2(v.z, t.x));
    out_support_and_tangent_vector[2] = packSnorm2x16(t.yz);
}