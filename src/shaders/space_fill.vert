#version 330 core

#ifndef ATOM_VEL
#define ATOM_VEL 1
#endif

#ifndef ATOM_COL
#define ATOM_COL 1
#endif

#ifndef ATOM_IDX
#define ATOM_IDX 1
#endif

layout (location = 0) in vec3  in_atom_pos;
layout (location = 1) in vec3  in_atom_prev_pos;
layout (location = 2) in float in_atom_rad;
layout (location = 3) in uint  in_atom_flags;
layout (location = 4) in vec4  in_atom_col;

out VS_GS {
    flat vec4 view_sphere;
#if ATOM_VEL
    flat vec3 view_velocity;
#endif
#if ATOM_COL
    flat vec4 color;
#endif
#if ATOM_IDX
    flat uint atom_idx;
#endif
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
    
    float u_radius_scale;
};

void main() {
	vec4  col = in_atom_col;
	vec3  pos = in_atom_pos;
    float rad = in_atom_rad * u_radius_scale;
    if (u_atom_mask > 0U && (in_atom_flags & u_atom_mask) == 0U) rad = 0;

    vec3 view_pos    = vec3(u_world_to_view * vec4(pos, 1.0));
	vec4 view_sphere = vec4(view_pos, rad);

	out_vert.view_sphere = view_sphere;
#if ATOM_VEL
	out_vert.view_velocity = vec3(u_world_to_view * vec4(in_atom_pos - in_atom_prev_pos, 0));
#endif
#if ATOM_COL
	out_vert.color = col;
#endif
#if ATOM_IDX
	out_vert.atom_idx = uint(gl_VertexID);
#endif
}