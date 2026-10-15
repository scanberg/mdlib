#version 330 core

// Particles derived from the atoms of another molecule, their source (md_gl_mol_set_derived): each is the weighted sum
// of its atoms, every atom taken in its periodic image nearest to the reference atom of the particle. Its velocity is the
// same sum of theirs. Drawn as one point per particle and captured with transform feedback (separate attributes) into
// the positions and velocities of the particles.

uniform samplerBuffer  u_buf_atom_pos;      // Positions of the atoms of the source (vec3)
uniform samplerBuffer  u_buf_atom_vel;      // Velocities of the atoms of the source (vec3)
uniform usamplerBuffer u_buf_derived_atom;  // The atoms of the particles, one after the other
uniform samplerBuffer  u_buf_derived_weight;// Their weights

uniform mat3 u_cell;        // The basis vectors of the unit cell, as columns
uniform mat3 u_inv_cell;
uniform vec3 u_pbc_mask;    // 1 along the periodic basis vectors, 0 along the others (all 0 without a cell)

layout (location = 0) in uvec3 in_derived;  // The offset of the atoms of the particle, their count, its reference atom

out vec3 out_position;
out vec3 out_velocity;

void main() {
    uint beg   = in_derived.x;
    uint count = in_derived.y;
    vec3 ref   = texelFetch(u_buf_atom_pos, int(in_derived.z)).xyz;

    // Relative to the reference, which keeps the precision of a particle far from the origin
    vec3  sum_d = vec3(0.0);
    vec3  sum_v = vec3(0.0);
    float sum_w = 0.0;
    for (uint i = 0u; i < count; ++i) {
        int   j = int(beg + i);
        int   a = int(texelFetch(u_buf_derived_atom, j).x);
        float w = texelFetch(u_buf_derived_weight, j).x;
        vec3  d = texelFetch(u_buf_atom_pos, a).xyz - ref;
        vec3  f = u_inv_cell * d;
        d -= u_cell * (floor(f + 0.5) * u_pbc_mask);
        sum_d += w * d;
        sum_v += w * texelFetch(u_buf_atom_vel, a).xyz;
        sum_w += w;
    }

    out_position = ref * sum_w + sum_d;
    out_velocity = sum_v;
}
