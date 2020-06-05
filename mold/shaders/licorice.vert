#version 330 core

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
    
    float u_radius;
};

layout(location = 0) in vec3 in_position;
layout(location = 1) in vec3 in_prev_position;
layout(location = 2) in uint in_atom_flags;
layout(location = 3) in vec4 in_color;

out Vertex {
    flat vec3 view_velocity;
    flat uint picking_idx;
    flat vec4 color;
} out_vert;

void main() {
    vec4 color = in_color;
    if (u_atom_mask > 0U && (in_atom_flags & u_atom_mask) == 0U) color = vec4(0);

    gl_Position = u_world_to_view * vec4(in_position, 1.0);
    out_vert.view_velocity = vec3(u_world_to_view * vec4(in_position - in_prev_position, 0.0));
    out_vert.picking_idx = uint(gl_VertexID);
    out_vert.color = color;
}