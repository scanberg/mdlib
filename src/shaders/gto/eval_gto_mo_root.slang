// eval_gto_mo_root.slang
// Simple root-argument variant of the MO evaluator.
// Buffer arguments are physical device pointers passed through push constants;
// the output image intentionally remains slot-bound for now.

#define OP_SET 0
#define OP_ADD 1
#define OP_SUB 2
#define OP_MAX 3
#define OP_MIN 4
#define OP_ACCUM_MASK 0x0000000Fu
#define OP_ABS 0x00000100u

struct PGTO {
    float coeff;
    float alpha;
    uint  ijkl;
};

struct RootArgs {
    float4x4 world_to_model;
    float4x4 index_to_world;
    float4   step;
    uint4    grid_dim;
    uint     num_cgtos;
    uint     num_rows;
    uint     mode;
    uint     operation;

    uint*   cgto_atom_idx;
    float*  cgto_r;
    uint2*  cgto_off_len;
    PGTO*   pgto;
    float4* atom_xyz;
    float*  coeffs;
};

[[vk::push_constant]]
ConstantBuffer<RootArgs> args;

[[vk::binding(16)]]
RWTexture3D<float> out_vol;

uint4 unpack_ijkl(uint packed) {
    return uint4(
        (packed      ) & 0xFFu,
        (packed >>  8) & 0xFFu,
        (packed >> 16) & 0xFFu,
        (packed >> 24) & 0xFFu
    );
}

float pow_upto4(float x, uint p) {
    float x2 = x * x;
    return (p == 0u) ? 1.0f :
           (p == 1u) ? x :
           (p == 2u) ? x2 :
           (p == 3u) ? x2 * x : x2 * x2;
}

float eval_cgto(uint cgto_idx, float3 coord) {
    float r = args.cgto_r[cgto_idx];
    if (r == 0.0f) {
        return 0.0f;
    }

    uint atom_idx = args.cgto_atom_idx[cgto_idx];
    float3 center = args.atom_xyz[atom_idx].xyz;
    float3 dc = coord - center;
    if (dot(dc, dc) > r * r) {
        return 0.0f;
    }

    uint2 off_len = args.cgto_off_len[cgto_idx];
    float value = 0.0f;
    for (uint i = 0u; i < off_len.y; ++i) {
        PGTO prim = args.pgto[off_len.x + i];
        uint4 ijkl = unpack_ijkl(prim.ijkl);
        float3 d = coord - center;
        float r2 = dot(d, d);
        value += prim.coeff *
            pow_upto4(d.x, ijkl.x) *
            pow_upto4(d.y, ijkl.y) *
            pow_upto4(d.z, ijkl.z) *
            exp(-prim.alpha * r2);
    }
    return value;
}

[numthreads(8, 8, 8)]
void main(uint3 dtid: SV_DispatchThreadID) {
    if (dtid.x >= args.grid_dim.x || dtid.y >= args.grid_dim.y || dtid.z >= args.grid_dim.z) {
        return;
    }

    float3 coord = mul(args.index_to_world, float4(float3(dtid), 1.0f)).xyz;
    float result = 0.0f;

    for (uint row = 0u; row < args.num_rows; ++row) {
        float psi = 0.0f;
        for (uint cgto_idx = 0u; cgto_idx < args.num_cgtos; ++cgto_idx) {
            float phi = eval_cgto(cgto_idx, coord);
            psi += args.coeffs[row * args.num_cgtos + cgto_idx] * phi;
        }
        result += (args.mode != 0u) ? psi * psi : psi;
    }

    if ((args.operation & OP_ABS) != 0u) {
        result = abs(result);
    }

    switch (args.operation & OP_ACCUM_MASK) {
        case OP_SET: break;
        case OP_ADD: result += out_vol[dtid]; break;
        case OP_SUB: result = out_vol[dtid] - result; break;
        case OP_MAX: result = max(result, out_vol[dtid]); break;
        case OP_MIN: result = min(result, out_vol[dtid]); break;
        default: break;
    }

    out_vol[dtid] = result;
}
