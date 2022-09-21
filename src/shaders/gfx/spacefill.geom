#version 460

#include "common.h"

#ifndef ORTHO
#define ORTHO 0
#endif

layout (binding = 0, std140) uniform UboBuffer {
    UniformData ubo;
};

layout (points) in;
layout (triangle_strip, max_vertices = 4) out;

layout(location = 0) in VS_GS {
    flat vec4 view_sphere;
    flat vec3 view_velocity;
    flat uint color;
    flat uint index;
} IN[];

layout(location = 0) out GS_FS {
    flat   vec4 view_sphere;
    flat   vec3 view_velocity;
    flat   uint color;
    flat   uint index;
    smooth vec3 view_coord;
#if ORTHO
    smooth vec2 uv;
#endif
} OUT;

#if ORTHO
void emit_vertex_ortho(vec2 uv) {
    vec3  pos = IN[0].view_sphere.xyz;
    float rad = IN[0].view_sphere.w;
    vec3  view_coord = vec3(pos + vec3(uv * rad, 0));

    OUT.view_sphere   = IN[0].view_sphere;
    OUT.view_velocity = IN[0].view_velocity;
    OUT.color         = IN[0].color;
    OUT.index         = IN[0].index;
    OUT.view_coord    = view_coord;
    OUT.uv            = uv;
    gl_Position       = ubo.view_to_clip * vec4(view_coord + vec3(0,0,rad), 1.0); 
    EmitVertex();
}
#else
void emit_vertex_persp(vec4 clip_coord) {
    vec4 view_coord = ubo.clip_to_view * clip_coord;

    OUT.view_sphere   = IN[0].view_sphere;
    OUT.view_velocity = IN[0].view_velocity;
    OUT.color         = IN[0].color;
    OUT.index         = IN[0].index;
    OUT.view_coord    = view_coord.xyz / view_coord.w;
    gl_Position       = clip_coord;
    EmitVertex();
}
#endif

// From Inigo Quilez!
void proj_sphere(in vec4 sphere, 
                 in vec2 fle,
                 out vec2 axis_a,
                 out vec2 axis_b,
                 out vec2 center) {
    vec3  o = sphere.xyz;
    float r2 = sphere.w*sphere.w;
    float z2 = o.z*o.z; 
    float l2 = dot(o,o);
    float c = -r2*(r2-l2)/((l2-z2)*(r2-z2));
    
    // axis
    axis_a = fle*sqrt(c/(r2-z2)) * vec2( o.x,o.y);
    axis_b = fle*sqrt(c/(r2-l2)) * vec2(-o.y,o.x);
    center = -fle*o.z*o.xy/(z2-r2);
}

void main() {
#if ORTHO
    emit_vertex_ortho(vec2(-1, -1));
    emit_vertex_ortho(vec2(+1, -1));
    emit_vertex_ortho(vec2(-1, +1));
    emit_vertex_ortho(vec2(+1, +1));
    EndPrimitive();
#else
    // Focal length
    float clip_z = -ubo.view_to_clip[2][2] - ubo.view_to_clip[3][2] / (IN[0].view_sphere.z + IN[0].view_sphere.w);
    vec2 fle = vec2(ubo.view_to_clip[0][0], ubo.view_to_clip[1][1]);
    vec2 axis_a;
    vec2 axis_b;
    vec2 center;
    proj_sphere(IN[0].view_sphere, fle, axis_a, axis_b, center);

    emit_vertex_persp(vec4(center -axis_a -axis_b, clip_z, 1));
    emit_vertex_persp(vec4(center +axis_a -axis_b, clip_z, 1));
    emit_vertex_persp(vec4(center -axis_a +axis_b, clip_z, 1));
    emit_vertex_persp(vec4(center +axis_a +axis_b, clip_z, 1));
    EndPrimitive();
#endif
}