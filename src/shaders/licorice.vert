#version 330 core

#ifndef ATOM_COL
#define ATOM_COL 1
#endif

#ifndef ATOM_VEL
#define ATOM_VEL 0
#endif

#ifndef ATOM_IDX
#define ATOM_IDX 1
#endif

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
layout(location = 1) in vec3 in_prev_position;
layout(location = 2) in uint in_atom_flags;
layout(location = 3) in vec4 in_color;

out Vertex {
#if ATOM_VEL
    flat vec3 view_velocity;
#endif
#if ATOM_IDX
    flat uint picking_idx;
#endif
#if ATOM_COL
    flat vec4 color;
#endif
    flat uint flags;
} out_vert;

void main() {
    gl_Position = u_world_to_view * vec4(in_position, 1.0);
    
    out_vert.flags = in_atom_flags;
#if ATOM_COL
    out_vert.color = in_color;
#endif
#if ATOM_VEL
    out_vert.view_velocity = vec3(u_world_to_view * vec4(in_position - in_prev_position, 0.0));
#endif
#if ATOM_IDX
    out_vert.picking_idx = uint(gl_VertexID);
#endif
}