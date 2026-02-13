#version 330 core
#extension GL_ARB_shading_language_packing : enable

uniform samplerBuffer u_buf_atom_pos;
uniform samplerBuffer u_buf_atom_vel;
uniform samplerBuffer u_buf_secondary_structure;

layout (location = 0) in uint in_bb_seg_idx;
layout (location = 1) in uint in_atom_offset;
layout (location = 2) in uint in_atom_ca;
layout (location = 3) in uint in_atom_c;
layout (location = 4) in uint in_atom_o;
layout (location = 5) in uint in_flags;

out vec3  out_position;
out uint  out_atom_idx;
out vec3  out_velocity;
out float out_segment_t;
out uint  out_secondary_structure_and_flags;
out uvec3 out_support_and_tangent_vector;

vec3 safe_normalize(vec3 v) {
    float d2 = dot(v, v);
    return d2 > 0.0 ? v * inversesqrt(d2) : vec3(0.0, 0.0, 1.0);
}

void main() {
    uint ca_idx = in_atom_offset + in_atom_ca;
    uint c_idx  = in_atom_offset + in_atom_c;
    uint o_idx  = in_atom_offset + in_atom_o;

    vec3 ca   = texelFetch(u_buf_atom_pos, int(ca_idx)).xyz;
    vec3 ca_v = texelFetch(u_buf_atom_vel, int(ca_idx)).xyz;
    vec3 c    = texelFetch(u_buf_atom_pos, int(c_idx)).xyz;
    vec3 o    = texelFetch(u_buf_atom_pos, int(o_idx)).xyz;
    vec3 ss   = texelFetch(u_buf_secondary_structure, int(in_bb_seg_idx)).xyz;

    vec3 support = safe_normalize(o - c);
    vec3 tangent = vec3(0.0);

    out_position = ca;
    out_atom_idx = ca_idx;
    out_velocity = ca_v;
    out_segment_t = 0.0;
    out_secondary_structure_and_flags = (in_flags << 24U) | (packUnorm4x8(vec4(ss, 0.0)) & 0x00FFFFFFU);
    out_support_and_tangent_vector[0] = packSnorm2x16(support.xy);
    out_support_and_tangent_vector[1] = packSnorm2x16(vec2(support.z, tangent.x));
    out_support_and_tangent_vector[2] = packSnorm2x16(tangent.yz);
}
