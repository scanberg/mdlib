#version 430 core

#ifndef GROUP_SIZE
#define GROUP_SIZE 16
#endif

layout (local_size_x = GROUP_SIZE) in;

struct AABB {
	float min_x, min_y, min_z;
	float max_x, max_y, max_z;
};

struct Position {
	float x, y, z;
};

struct ResidueRange {
	uint beg;
	uint end;
};

layout (std430, binding = 0) writeonly buffer aabb_buffer {
	uint out_aabb[];
};

layout (std430, binding = 1) readonly buffer position_buffer {
	Position in_position[];
};

layout (std430, binding = 2) readonly buffer radius_buffer {
	float in_radius[];
};

layout (std430, binding = 3) readonly buffer residue_buffer {
	ResidueRange in_range[];
};


void main() {
	uint g_idx = gl_LocalInvocationIndex;
	uint l_idx = gl_LocalInvocationID.x;

	ResidueRange range = in_range[g_idx];
}