#version 460

#include "common.h"
#include "culling.glsl"

#ifndef ORTHO
#define ORTHO 0
#endif

layout (binding = 0, std140) uniform UboBuffer {
    UniformData ubo;
};

layout(binding = 0) uniform sampler2D in_depth_pyramid_tex;

layout (points) in;
layout (triangle_strip, max_vertices = 4) out;

layout(location = 0) in VS_GS {
    flat vec4 view_sphere;
    flat uint index;
} IN[];

layout(location = 0) out GS_FS {
    flat   vec4 view_sphere;
    smooth vec3 view_coord;
    flat   uint index;
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

// 2D Polyhedral Bounds of a Clipped, Perspective-Projected 3D Sphere. Michael Mara, Morgan McGuire. 2013
// There is an assumption here that the projection matrix is a symmetric projection matrix, hence only P[0][0] and P[1][1] is used.
void project_sphere(vec4 S, float P00, float P11, out vec4 aabb)
{
    P00 = -P00;
    P11 = -P11;
	vec2 cx = -S.xz;
	vec2 vx = vec2(sqrt(dot(cx, cx) - S.w * S.w), S.w);
	vec2 minx = mat2(vx.x, vx.y, -vx.y, vx.x) * cx;
	vec2 maxx = mat2(vx.x, -vx.y, vx.y, vx.x) * cx;

	vec2 cy = -S.yz;
	vec2 vy = vec2(sqrt(dot(cy, cy) - S.w * S.w), S.w);
	vec2 miny = mat2(vy.x, vy.y, -vy.y, vy.x) * cy;
	vec2 maxy = mat2(vy.x, -vy.y, vy.y, vy.x) * cy;

	aabb = vec4(minx.x / minx.y * P00, miny.x / miny.y * P11, maxx.x / maxx.y * P00, maxy.x / maxy.y * P11);
}

#define AABB 1

void main() {
#if ORTHO
    emit_vertex_ortho(vec2(-1, -1));
    emit_vertex_ortho(vec2(+1, -1));
    emit_vertex_ortho(vec2(-1, +1));
    emit_vertex_ortho(vec2(+1, +1));
    EndPrimitive();
#else
    mat4 P = ubo.view_to_clip;

    float clip_z = -P[2][2] - P[3][2] / (IN[0].view_sphere.z + IN[0].view_sphere.w);
    #if AABB
        vec4 aabb;
        project_sphere(IN[0].view_sphere, P[0][0], P[1][1], aabb);

        bool visible = true;

//#if ENABLE_CULLING
		if (ubo.late == 1 && ubo.depth_culling == 1) {
		    vec3 aabb_uv_min = vec3(aabb.xy, clip_z) * 0.5 + 0.5;
		    vec3 aabb_uv_max = vec3(aabb.zw, clip_z) * 0.5 + 0.5;

            visible = cull_hiz(aabb_uv_min, aabb_uv_max, in_depth_pyramid_tex, ubo.depth_pyramid_width, ubo.depth_pyramid_height);
		}
//#endif

        if (visible) {
            emit_vertex_persp(vec4(aabb.xy, clip_z, 1));
            emit_vertex_persp(vec4(aabb.zy, clip_z, 1));
            emit_vertex_persp(vec4(aabb.xw, clip_z, 1));
            emit_vertex_persp(vec4(aabb.zw, clip_z, 1));
        }

        EndPrimitive();
    #else
        // Focal length
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
#endif
}