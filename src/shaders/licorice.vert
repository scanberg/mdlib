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
    uint _pad0;
    uint _pad1;
    uint _pad2;
    float u_radius;
};

layout(location = 0) in vec3 in_position;
layout(location = 1) in vec3 in_velocity;
layout(location = 2) in uint in_atom_flags;
layout(location = 3) in vec4 in_color;

out Vertex {
    flat vec3 view_vel;
    flat uint picking_idx;
    flat vec4 color;
    flat uint flags;
} out_vert;

void main() {
    out_vert.flags = in_atom_flags;
    out_vert.color = in_color;
    out_vert.view_vel = vec3(u_world_to_view * vec4(in_velocity, 0.0));
    out_vert.picking_idx = uint(gl_VertexID);
    gl_Position = u_world_to_view * vec4(in_position, 1.0);
}