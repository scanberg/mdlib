#version 330 core

layout (location = 0) in vec3  in_position;
layout (location = 1) in uint  in_atom_index;
layout (location = 2) in vec3  in_velocity;
layout (location = 3) in float in_segment_t;
layout (location = 4) in vec3  in_secondary_structure;
layout (location = 5) in uint  in_flags;
layout (location = 6) in vec3  in_support_vector;

out Vertex {
    vec3  position;
    uint  atom_index;
    vec3  velocity;
    float segment_t;
    vec3  secondary_structure;
    uint  flags;
    vec3  support_vector;
} out_vert;

void main() {
    out_vert.position 		 	 = in_position;
    out_vert.atom_index 		 = in_atom_index;
    out_vert.velocity 			 = in_velocity;
    out_vert.segment_t 			 = in_segment_t;
    out_vert.secondary_structure = in_secondary_structure;
    out_vert.flags 				 = in_flags;
    out_vert.support_vector      = in_support_vector;
} 