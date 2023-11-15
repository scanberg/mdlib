#version 150 core

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
    uint _pad0;
    float u_radius;
};

layout (lines) in;
layout (triangle_strip, max_vertices = 24) out;

in Vertex {
    flat vec3 view_vel;
    flat uint picking_idx;
    flat vec4 color;
    flat uint flags;
} in_vert[];

out Fragment {
    flat uint picking_idx[2];
    flat vec4 color[2];
    flat vec4 capsule_center_radius;
    flat vec4 capsule_axis_length;
    smooth vec3 view_vel;
    smooth vec3 view_pos;
} out_frag;

vec4 view_vertices[8];
vec4 proj_vertices[8];

void emit_vertex(int i) {
    out_frag.view_vel = in_vert[i/4].view_vel;
    out_frag.view_pos = view_vertices[i].xyz;
    gl_Position = proj_vertices[i];
    EmitVertex();
}

void emit(int a, int b, int c, int d) {
    emit_vertex(a);
    emit_vertex(b);
    emit_vertex(c);
    emit_vertex(d);
    EndPrimitive(); 
}

vec3 get_ortho_vec(vec3 v, vec3 A, vec3 B) {
    if(abs(1-dot(v,A))>0.001){
        return normalize(cross(v,A));
    }else{
        return normalize(cross(v,B));
    }
}

void main() {
    if ((in_vert[0].flags & u_atom_mask) != u_atom_mask ||
        (in_vert[1].flags & u_atom_mask) != u_atom_mask) {
        return;
    }

    if (in_vert[0].color.a == 0 || in_vert[1].color.a == 0) {
        return;
    }

    // Compute orientation vectors for the two connecting faces:
    vec3 p0 = gl_in[0].gl_Position.xyz;
    vec3 p1 = gl_in[1].gl_Position.xyz;
    float r = u_radius;
    float l = distance(p0, p1);
    vec3 a = (p1 - p0) / l;
    vec3 c = (p0 + p1) * 0.5;

    out_frag.color[0] = in_vert[0].color;
    out_frag.color[1] = in_vert[1].color;

    out_frag.picking_idx[0] = u_bond_base_index + uint(gl_PrimitiveIDIn);//in_vert[0].picking_idx;
    out_frag.picking_idx[1] = u_bond_base_index + uint(gl_PrimitiveIDIn);//in_vert[1].picking_idx;

    out_frag.capsule_center_radius = vec4(c, r);
    out_frag.capsule_axis_length = vec4(a, l);

    // Extend end points to properly fit the sphere caps
    p0 -= a * r;
    p1 += a * r;

    const vec3 B = vec3(0,0,1);
    const vec3 A = vec3(1,0,0);
    vec3 o = get_ortho_vec(a,A,B);

    // Declare scratch variables for basis vectors:
    vec3 i,j,k;

    // Compute vertices of prismoid:
    j = a; i = o; k = normalize(cross(i, j)); i = cross(k, j); i *= r; k *= r;
    view_vertices[0] = vec4(p0 + i + k, 1);
    view_vertices[1] = vec4(p0 + i - k, 1);
    view_vertices[2] = vec4(p0 - i - k, 1);
    view_vertices[3] = vec4(p0 - i + k, 1);
    view_vertices[4] = vec4(p1 + i + k, 1);
    view_vertices[5] = vec4(p1 + i - k, 1);
    view_vertices[6] = vec4(p1 - i - k, 1);
    view_vertices[7] = vec4(p1 - i + k, 1);

    proj_vertices[0] = u_view_to_clip * view_vertices[0];
    proj_vertices[1] = u_view_to_clip * view_vertices[1];
    proj_vertices[2] = u_view_to_clip * view_vertices[2];
    proj_vertices[3] = u_view_to_clip * view_vertices[3];
    proj_vertices[4] = u_view_to_clip * view_vertices[4];
    proj_vertices[5] = u_view_to_clip * view_vertices[5];
    proj_vertices[6] = u_view_to_clip * view_vertices[6];
    proj_vertices[7] = u_view_to_clip * view_vertices[7];

    // Emit the six faces of the prismoid:
    emit(0,1,3,2); emit(5,4,6,7);
    emit(4,5,0,1); emit(3,2,7,6);
    emit(0,3,4,7); emit(2,1,6,5);
}