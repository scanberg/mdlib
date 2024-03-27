#version 330 core

#ifndef ORTHO
#define ORTHO 0
#endif

layout (points) in;
layout (triangle_strip, max_vertices = 4) out;

in VS_GS {
    flat vec4 view_sphere;
    flat vec3 view_velocity;
    flat vec4 color;
    flat uint atom_idx;
} in_vert[];

out GS_FS {
    smooth vec3 view_coord;
    flat vec4 view_sphere;
    flat vec3 view_velocity;
    flat vec4 color;
    flat uint atom_idx;
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
    uint _pad0;
    uint _pad1;
    uint _pad2;
    float u_radius_scale;
};

void emit_vertex_ortho(vec2 uv) {
    vec3  pos = in_vert[0].view_sphere.xyz;
    float rad = in_vert[0].view_sphere.w;
    out_frag.uv = uv;
    out_frag.view_coord = vec3(pos + vec3(uv * rad, 0));
    gl_Position = u_view_to_clip * vec4(out_frag.view_coord + vec3(0,0,rad), 1.0);
    EmitVertex();
}

void emit_vertex_persp(vec4 clip_coord) {
    vec4 view_coord = u_clip_to_view * clip_coord;
    out_frag.view_coord = view_coord.xyz / view_coord.w;
    gl_Position = clip_coord;
    EmitVertex();
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

void main()
{
    if (in_vert[0].view_sphere.w == 0) {
        EndPrimitive();
        return;
    }
    if (in_vert[0].color.a == 0) {
        EndPrimitive();
        return;
    }

    out_frag.view_sphere = in_vert[0].view_sphere;
    out_frag.view_velocity = in_vert[0].view_velocity;
    out_frag.color = in_vert[0].color;
    out_frag.atom_idx = in_vert[0].atom_idx;

#if ORTHO
    emit_vertex_ortho(vec2(-1, -1));
    emit_vertex_ortho(vec2(+1, -1));
    emit_vertex_ortho(vec2(-1, +1));
    emit_vertex_ortho(vec2(+1, +1));
    EndPrimitive();
#else
    // Focal length
    float clip_z = -u_view_to_clip[2][2] - u_view_to_clip[3][2] / (in_vert[0].view_sphere.z + in_vert[0].view_sphere.w);
    vec4 aabb = project_sphere(in_vert[0].view_sphere, u_view_to_clip);

    emit_vertex_persp(vec4(aabb.x, aabb.z, clip_z, 1));
    emit_vertex_persp(vec4(aabb.y, aabb.z, clip_z, 1));
    emit_vertex_persp(vec4(aabb.x, aabb.w, clip_z, 1));
    emit_vertex_persp(vec4(aabb.y, aabb.w, clip_z, 1));

    EndPrimitive();
#endif
}