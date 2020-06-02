#version 330 core

layout (location = 0) in uvec2 in_atom_range;

out vec3 out_aabb_min;
out vec3 out_aabb_max;

uniform samplerBuffer u_atom_pos_buf;
uniform samplerBuffer u_atom_rad_buf;

void main() {
	int beg = int(in_atom_range.x);
	int end = int(in_atom_range.y);

	vec3  pos = texelFetch(u_atom_pos_buf, beg).xyz;
	float rad = texelFetch(u_atom_rad_buf, beg).r; 
	vec3 min_box = pos - rad;
	vec3 max_box = pos + rad;

	for (int i = beg + 1; i < end; i++) {
		pos = texelFetch(u_atom_pos_buf, i).xyz;
		rad = texelFetch(u_atom_rad_buf, i).r;
		min_box = min(min_box, pos - rad);
		max_box = max(max_box, pos + rad);
	}

	out_aabb_min = min_box;
	out_aabb_max = max_box;
}