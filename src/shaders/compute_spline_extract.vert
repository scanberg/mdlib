#version 330 core
#extension GL_ARB_shading_language_packing : enable

// The control points of a backbone from the atoms of their residues: the CA, and the carbonyl (O - C) as the fallback
// orientation (compute_spline_orient.vert), or with FLAG_FIXED_FACE the direction the flat faces of the profile face.
// The secondary structure (coil, helix, sheet) is in the xyz of u_buf_secondary_structure; its w is the radius of a
// tube (md_gl_mol_set_tubes) as a fraction of TUBE_RADIUS_MAX, 0 for a protein, and its y the half length of the flat
// sides of a tube. A tube is the backbone of a molecule of particles, each control point a particle (in ca).

// Must match md_gl.c, compute_spline_subdivide.vert and backbone.vert
#define TUBE_RADIUS_MAX 8.0

uniform samplerBuffer u_buf_atom_pos;
uniform samplerBuffer u_buf_atom_vel;
uniform samplerBuffer u_buf_secondary_structure;

layout (location = 0) in uint in_bb_seg_idx;
layout (location = 1) in uint in_atom_ca;
layout (location = 2) in uint in_atom_c;
layout (location = 3) in uint in_atom_o;
layout (location = 4) in uint in_flags;

out vec3  out_position;
out uint  out_atom_idx;
out vec3  out_velocity;
out float out_radius;
out uint  out_secondary_structure_and_flags;
out uvec3 out_support_and_tangent_vector;

vec3 safe_normalize(vec3 v) {
    float d2 = dot(v, v);
    return d2 > 0.0 ? v * inversesqrt(d2) : vec3(0.0, 0.0, 1.0);
}

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
    uint ca_idx = in_atom_ca;
    uint c_idx  = in_atom_c;
    uint o_idx  = in_atom_o;

    vec3 ca   = texelFetch(u_buf_atom_pos, int(ca_idx)).xyz;
    vec3 ca_v = texelFetch(u_buf_atom_vel, int(ca_idx)).xyz;
    vec3 c    = texelFetch(u_buf_atom_pos, int(c_idx)).xyz;
    vec3 o    = texelFetch(u_buf_atom_pos, int(o_idx)).xyz;
    vec4 ss   = texelFetch(u_buf_secondary_structure, int(in_bb_seg_idx));

    vec3 support = safe_normalize(o - c);
    vec3 tangent = vec3(0.0);

    out_position = ca;
    out_atom_idx = ca_idx;
    out_velocity = ca_v;
    out_radius = ss.w * TUBE_RADIUS_MAX;
    out_secondary_structure_and_flags = ((in_flags & 0xFFU) << 24U) | (packUnorm4x8(vec4(ss.xyz, 0.0)) & 0x00FFFFFFU);
    out_support_and_tangent_vector[0] = packSnorm2x16(support.xy);
    out_support_and_tangent_vector[1] = packSnorm2x16(vec2(support.z, tangent.x));
    out_support_and_tangent_vector[2] = packSnorm2x16(tangent.yz);
}
