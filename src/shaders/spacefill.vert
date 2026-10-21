#version 330 core

#ifndef ORTHO
#define ORTHO 0
#endif

// SPACE FILL. Every atom is one quad of 4 vertices that bounds the image of its sphere, drawn by vertex pulling without
// a geometry shader: atom a = quad_primitive() (QUADS in md_gl.c), corner c = gl_VertexID & 3 (bit 0: right, bit 1:
// top), from the static index buffer of the quads, two triangles (0, 1, 2), (2, 1, 3) each. The quad lies in front of
// the sphere, so the depth spacefill.frag writes only increases (layout (depth_greater)), and holds the points of the
// image plane whose rays the fragment shader intersects with it:
//   - perspective: the bounds of the projected sphere (Mara and McGuire 2013) on the plane of its nearest point;
//   - orthographic: the square of its silhouette, moved to its nearest point.

out VS_FS {
    smooth vec3 view_coord;
    flat vec4 view_sphere;
    flat vec3 view_velocity;
    flat vec4 color;
    flat uint atom_idx;
    flat uint flags;
    smooth vec2 uv;
} out_frag;

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

    float u_radius_scale;
};

uniform samplerBuffer  u_buf_atom_pos;      // RGB32F
uniform samplerBuffer  u_buf_atom_vel;      // RGB32F
uniform samplerBuffer  u_buf_atom_rad;      // R32F
uniform usamplerBuffer u_buf_atom_flags;    // R8UI
uniform samplerBuffer  u_buf_atom_color;    // RGBA8

// The instances of the draw (md_gl.c): the copies, the blocks of quads, the lists
#pragma INSTANCES

void cull() {
    // All corners of the quad coincide: no area, no fragments
    gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
}

// 2D Polyhedral Bounds of a Clipped, Perspective-Projected 3D Sphere. Michael Mara, Morgan McGuire. 2013
// Should work for asymmetric projection matrices
vec4 project_sphere(vec4 S, mat4 P) {
    vec2 cx = vec2(S.x, -S.z);
    vec2 vx = vec2(sqrt(dot(cx, cx) - S.w * S.w), S.w);
    vec2 minx = mat2(vx.x, -vx.y, vx.y, vx.x) * cx;
    vec2 maxx = mat2(vx.x, vx.y, -vx.y, vx.x) * cx;

    vec2 cy = vec2(S.y, -S.z);
    vec2 vy = vec2(sqrt(dot(cy, cy) - S.w * S.w), S.w);
    vec2 miny = mat2(vy.x, -vy.y, vy.y, vy.x) * cy;
    vec2 maxy = mat2(vy.x, vy.y, -vy.y, vy.x) * cy;

    return vec4(
        minx.x / minx.y * P[0][0] - P[2][0],
        maxx.x / maxx.y * P[0][0] - P[2][0],
        miny.x / miny.y * P[1][1] - P[2][1],
        maxy.x / maxy.y * P[1][1] - P[2][1]);
}

void main() {
    uint atom = quad_primitive();
    if (atom == QUAD_NONE) {
        cull();
        return;
    }
    int  a      = int(atom);
    uint flags  = texelFetch(u_buf_atom_flags, a).x;
    vec4 color  = texelFetch(u_buf_atom_color, a);
    float rad   = texelFetch(u_buf_atom_rad, a).x * u_radius_scale;
    // Not of the atom mask, or not shown
    if ((flags & u_atom_mask) != u_atom_mask || rad == 0.0 || color.a == 0.0) {
        cull();
        return;
    }

    vec3 atom_pos = texelFetch(u_buf_atom_pos, a).xyz;
    vec3 atom_vel = texelFetch(u_buf_atom_vel, a).xyz;
    vec3 pos      = instance_point(atom_pos);
    vec4 view_sphere = vec4(vec3(u_world_to_view * vec4(pos, 1.0)), rad);

    out_frag.view_sphere   = view_sphere;
    out_frag.view_velocity = vec3(u_world_to_view * vec4(instance_velocity(atom_pos, atom_vel), 0.0));
    out_frag.color         = color;
    out_frag.atom_idx      = u_atom_base_index + atom;
    out_frag.flags         = flags;

    uint corner = uint(gl_VertexID) & 3u;
    vec2 uv = vec2((corner & 1u) != 0u ? 1.0 : -1.0, (corner & 2u) != 0u ? 1.0 : -1.0);
    out_frag.uv = uv;

#if ORTHO
    out_frag.view_coord = view_sphere.xyz + vec3(uv * rad, 0.0);
    gl_Position = u_view_to_clip * vec4(out_frag.view_coord + vec3(0.0, 0.0, rad), 1.0);
#else
    // On the plane of the sphere's nearest point
    float clip_z = -u_view_to_clip[2][2] - u_view_to_clip[3][2] / (view_sphere.z + view_sphere.w);
    vec4  aabb   = project_sphere(view_sphere, u_view_to_clip);
    vec4  clip   = vec4(uv.x < 0.0 ? aabb.x : aabb.y, uv.y < 0.0 ? aabb.z : aabb.w, clip_z, 1.0);
    vec4  view   = u_clip_to_view * clip;
    out_frag.view_coord = view.xyz / view.w;
    gl_Position = clip;
#endif
}
