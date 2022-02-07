#version 330 core

uniform vec3 u_pbc_ext;

layout (location = 0) in vec3 in_position;
layout (location = 1) in vec3 in_position_prev;

out vec3 out_velocity;

vec3 de_periodize(vec3 pos, vec3 ref_pos, vec3 ext) {
    vec3 delta = pos - ref_pos;
    vec3 signed_mask = sign(delta) * step(ext * 0.5, abs(delta));
    return pos - ext * signed_mask;
}

void main() {
    vec3 cur = in_position;
    vec3 old = de_periodize(in_position_prev, in_position, u_pbc_ext);
    out_velocity = cur - old;
}