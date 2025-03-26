#version 330 core

uniform samplerBuffer u_buf_atom_pos;
uniform samplerBuffer u_buf_atom_vel;
uniform samplerBuffer u_buf_secondary_structure;

layout (location = 0) in uint in_bb_seg_idx;
layout (location = 1) in uint in_atom_offset;
layout (location = 2) in uint in_atom_ca;
layout (location = 3) in uint in_atom_c;
layout (location = 4) in uint in_atom_o;
layout (location = 5) in uint in_flags;

out Vertex {
    vec3  position;
    uint  atom_index;
    vec3  velocity;
    float segment_t;
    vec3  secondary_structure;
    uint  flags;
    vec3  support_vector;
} out_vert;

void main() {
    uint ca_idx = in_atom_offset + in_atom_ca;
    uint  c_idx = in_atom_offset + in_atom_c;
    uint  o_idx = in_atom_offset + in_atom_o;

    vec3 ca   = texelFetch(u_buf_atom_pos, int(ca_idx)).xyz;
    vec3 ca_v = texelFetch(u_buf_atom_vel, int(ca_idx)).xyz;
    vec3 c    = texelFetch(u_buf_atom_pos, int(c_idx)).xyz;
    vec3 o    = texelFetch(u_buf_atom_pos, int(o_idx)).xyz;
    vec3 ss   = texelFetch(u_buf_secondary_structure, int(in_bb_seg_idx)).xyz;
    uint f    = in_flags;

    vec3 v = normalize(o - c);

    out_vert.position 		 	 = ca;
    out_vert.atom_index 		 = ca_idx;
    out_vert.velocity 			 = ca_v;
    out_vert.segment_t 			 = 0;
    out_vert.secondary_structure = ss;
    out_vert.flags 				 = f;
    out_vert.support_vector      = v;
} 