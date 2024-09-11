#version 430 core

struct GTO {
    vec3  coord;
    float coeff;
    float alpha;
    float cutoff;
    uint  ijkl_packed;
};

layout (local_size_x = 8, local_size_y = 8, local_size_z = 8) in;
layout(rgba16f, binding = 0) uniform image3D out_vol;

layout(location = 0) uniform GTO  in_gtos[];
layout(location = 1) uniform uint in_num_gtos;

layout(location = 2) uniform vec3 in_offset;
layout(location = 3) uniform vec3 in_step_x;
layout(location = 4) uniform vec3 in_step_y;
layout(location = 5) uniform vec3 in_step_z;

shared uint num_gtos = 0;
shared GTO  gtos[256];

float fast_pow(float x, uint i) {
    return pow(x, float(i));
}

void main() {
    // Step 1: Load and test gtos against current block
    vec3 aabb_min;
    vec3 aabb_max;

    for (uint i = 0; i < in_num_gtos; ++i) {
        float cutoff = in_gtos[i].cutoff;
        vec3 coord   = in_gtos[i].coord;

        vec3 clamped = clamp(coord, aabb_min, aabb_max);
        vec3 delta = clamped - coord;
        float d2 = dot(delta, delta);

        if (d2 < cutoff * cutoff) {
            uint idx = atomicAdd(num_gtos, 1);
            out_gtos[idx] = in_gtos[i];
        }
    }
    
    // Step 2: Evaluate GTOs for current coord
    float psi = 0;

    vec3 coord = vec3(
        in_origin.x + ix * in_step_x.x + iy * in_step_y.x + iz * in_step_z.x,
        in_origin.y + ix * in_step_x.y + iy * in_step_y.y + iz * in_step_z.y,
        in_origin.z + ix * in_step_x.z + iy * in_step_y.z + iz * in_step_z.z
    );

    for (uint i = 0; i < num_gtos; ++i) {
        float px    = gtos[i].x;
        float py    = gtos[i].y;
        float pz    = gtos[i].z;
        float alpha = gtos[i].alpha;
        float coeff = gtos[i].coeff;
        uint  pi    = gtos[i].packed_ijkl & 0xFFu;
        uint  pj    = gtos[i].packed_ijkl >> 8u  & 0xFFu;
        uint  pk    = gtos[i].packed_ijkl >> 16u & 0xFFu;

        float dx = coord.x - px;
        float dy = coord.y - py;
        float dz = coord.z - pz;
        float d2 = dx * dx + dy * dy + dz * dz;
        float fx = fast_pow(dx, pi);
        float fy = fast_pow(dy, pj);
        float fz = fast_pow(dz, pk);
        float exp_term = exp(-alpha * d2);
        float powxyz = fx * fy * fz;
        float prod = coeff * powxyz * exp_term;
        psi += prod;
    }
    
    ivec3 vol_idx = gl_GlobalInvocationID.xyz;
    ivec3 vol_size = imageSize(out_vol);
    if (vol_idx.x < vol_size.x &&
        vol_idx.y < vol_size.y &&
        vol_idx.z < vol_size.z) {

        imageStore(out_vol, vol_idx, psi);
    }
}