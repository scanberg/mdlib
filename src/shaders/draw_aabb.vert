#version 330 core

layout(location = 0) in vec3 in_aabb_min;
layout(location = 1) in vec3 in_aabb_max;

out Vertex {
    flat vec3 aabb_min;
    flat vec3 aabb_ext;
    flat uint id;
} out_vert;

void main() {
    out_vert.aabb_min = in_aabb_min;
    out_vert.aabb_ext = in_aabb_max - in_aabb_min;
    out_vert.id = uint(gl_VertexID);
}