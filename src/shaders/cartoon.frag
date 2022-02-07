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
    float u_width_scale;
    float u_height_scale;
};


in Fragment {
    smooth vec3 view_coord;
    smooth vec3 view_velocity;
    smooth vec4 color;
    smooth vec3 view_normal;
    flat   uint picking_idx;
} in_frag;

#pragma EXTRA_SRC

void main() {
    vec4 color = in_frag.color;
    vec3 view_normal = normalize(in_frag.view_normal);
    vec3 view_coord  = in_frag.view_coord;
    vec3 view_velocity = in_frag.view_velocity;
    uint atom_index = in_frag.picking_idx;

    write_fragment(view_coord, view_velocity, view_normal, color, atom_index);
}