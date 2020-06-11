#version 330 core

uniform samplerBuffer u_atom_color_buffer;
uniform usamplerBuffer u_atom_flags_buffer;

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
    float u_width_scale;
    float u_height_scale;
};

layout (location = 0) in vec3  in_control_point;
layout (location = 1) in uint  in_atom_index;
layout (location = 2) in vec3  in_velocity;
layout (location = 3) in float in_segment_t;
layout (location = 4) in vec3  in_secondary_structure;
layout (location = 5) in uint  in_flags;
layout (location = 6) in vec3  in_support_vector;
layout (location = 7) in vec3  in_support_tangent;

out Vertex {
    vec4  control_point;
    vec4  support_vector;
    vec4  support_tangent;
    vec4  velocity;
    float segment_t;
    vec3  secondary_structure;
    uint  flags;
    vec4  color;
    uint  picking_idx;
} out_vert;
void main() {
    uint atom_flags = texelFetch(u_atom_flags_buffer, int(in_atom_index)).x;
    vec4 atom_color = texelFetch(u_atom_color_buffer, int(in_atom_index));
    if (u_atom_mask > 0U && (atom_flags & u_atom_mask) == 0U) atom_color = vec4(0);

    out_vert.control_point       = vec4(in_control_point, 1);
    out_vert.support_vector      = vec4(in_support_vector, 0);
    out_vert.support_tangent     = vec4(in_support_tangent, 0);
    out_vert.velocity            = vec4(in_velocity, 0);
    out_vert.segment_t           = in_segment_t;
    out_vert.secondary_structure = in_secondary_structure;
    out_vert.flags               = in_flags;
    out_vert.color               = atom_color;
    out_vert.picking_idx         = in_atom_index;
}