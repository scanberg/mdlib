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
    uint u_atom_base_index;
    uint u_bond_base_index;
    uint _pad;

    vec4  u_scale;
    uvec4 u_res;
    uvec4 u_rings;
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
    write_fragment(in_frag.view_coord, in_frag.view_velocity, normalize(in_frag.view_normal), in_frag.color, in_frag.picking_idx);
}
