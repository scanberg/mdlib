#version 330 core

layout (location = 0) in vec3  in_atom_pos;
layout (location = 1) in vec3  in_atom_vel;
layout (location = 2) in float in_atom_rad;
layout (location = 3) in uint  in_atom_flags;
layout (location = 4) in vec4  in_atom_col;

out VS_GS {
    flat vec4 view_sphere;
    flat vec3 view_velocity;
    flat vec4 color;
    flat uint atom_idx;
} out_vert;

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
    float u_radius_scale;
};

void main() {
	vec4  col = in_atom_col;
	vec3  pos = in_atom_pos;
    float rad = in_atom_rad * u_radius_scale;
    if ((in_atom_flags & u_atom_mask) != u_atom_mask) rad = 0;

    vec3 view_pos    = vec3(u_world_to_view * vec4(pos, 1.0));
	vec4 view_sphere = vec4(view_pos, rad);

	out_vert.view_sphere = view_sphere;
	out_vert.view_velocity = vec3(u_world_to_view * vec4(in_atom_vel, 0));
	out_vert.color = col;
	out_vert.atom_idx = uint(gl_VertexID);
}