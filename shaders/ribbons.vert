#version 330 core

uniform samplerBuffer u_atom_color_buffer;

layout (location = 0) in vec3 in_control_point;
layout (location = 1) in uint in_atom_index;
layout (location = 2) in vec3 in_velocity;
layout (location = 3) in uint in_flags;
layout (location = 4) in vec3 in_support_vector;
layout (location = 5) in vec3 in_support_tangent;

out Vertex {
    vec4 control_point;
    vec4 support_vector;
    vec4 support_tangent;
    vec4 velocity;
    uint flags;
    vec4 color;
    uint picking_idx;
} out_vert;

void main() {
    out_vert.control_point   = vec4(in_control_point, 1);
    out_vert.support_vector  = vec4(in_support_vector, 0);
    out_vert.support_tangent = vec4(in_support_tangent, 0);
    out_vert.velocity        = vec4(in_velocity, 0);
    out_vert.flags           = in_flags;
    out_vert.color           = texelFetch(u_atom_color_buffer, int(in_atom_index));
    out_vert.picking_idx     = in_atom_index;
}