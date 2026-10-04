// gto_common.slang
// Helpers shared by the shell-based GTO density kernels
// (eval_gto_density_tiled.slang and eval_gto_density_gemm.slang).
//
// Basis data on the device, per contracted shell (see md_gto.c, md_gto_gpu_basis_create):
//   shell_atom_idx[s]  : parent atom
//   shell_r[s]         : screening radius = max radius of influence over all components/primitives
//   shell_info[s]      : (prim_offset, num_prims, l, ao_offset)
//   shell_prim[p]      : (coeff, alpha * log2(e))  -> radial(r2) = sum coeff * exp2(-alpha' * r2)
// A Cartesian AO (shell s, component c) evaluates to
//   phi = norm(i,j,k) * x^i y^j z^k * radial(r2)
// with (i,j,k) in the canonical order documented in md_gto.h and
//   norm(i,j,k) = 1 / sqrt((2i-1)!! (2j-1)!! (2k-1)!!).

#ifndef GTO_COMMON_SLANG
#define GTO_COMMON_SLANG

#define OP_SET 0
#define OP_ADD 1
#define OP_SUB 2
#define OP_MAX 3
#define OP_MIN 4
#define OP_ACCUM_MASK 0x0000000Fu
#define OP_ABS 0x00000100u

static const float GTO_DF_ODD[5] = { 1.0f, 1.0f, 3.0f, 15.0f, 105.0f }; // (2n-1)!! for n = 0..4

uint gto_num_cart(uint l) {
    return ((l + 1u) * (l + 2u)) / 2u;
}

// (i,j,k) of component c of a shell with angular momentum l, packed 4 bits each.
uint gto_cart_ijk_packed(uint l, uint c) {
    uint idx = 0u;
    for (int i = int(l); i >= 0; --i) {
        for (int j = int(l) - i; j >= 0; --j) {
            if (idx == c) {
                uint k = l - uint(i) - uint(j);
                return uint(i) | (uint(j) << 4) | (k << 8);
            }
            idx++;
        }
    }
    return 0u;
}

float gto_cart_norm(uint ijk) {
    uint i = ijk & 0xFu;
    uint j = (ijk >> 4) & 0xFu;
    uint k = (ijk >> 8) & 0xFu;
    return rsqrt(GTO_DF_ODD[min(i, 4u)] * GTO_DF_ODD[min(j, 4u)] * GTO_DF_ODD[min(k, 4u)]);
}

// x^n for n <= 4 given x and x^2.
float gto_pow_le4(float x, float x2, uint n) {
    float r = 1.0f;
    r = (n == 1u) ? x       : r;
    r = (n == 2u) ? x2      : r;
    r = (n == 3u) ? x2 * x  : r;
    r = (n == 4u) ? x2 * x2 : r;
    return r;
}

float gto_monomial(float3 d, uint ijk) {
    float fx = gto_pow_le4(d.x, d.x * d.x,  ijk        & 0xFu);
    float fy = gto_pow_le4(d.y, d.y * d.y, (ijk >> 4)  & 0xFu);
    float fz = gto_pow_le4(d.z, d.z * d.z, (ijk >> 8)  & 0xFu);
    return fx * fy * fz;
}

// Packed upper-triangular index for a symmetric N x N matrix.
uint gto_tri_index(uint i, uint j, uint N) {
    uint row = min(i, j);
    uint col = max(i, j);
    return row * (2u * N - row + 1u) / 2u + (col - row);
}

float gto_apply_op(float value, float prev, uint operation) {
    if ((operation & OP_ABS) != 0u) value = abs(value);
    switch (operation & OP_ACCUM_MASK) {
        case OP_ADD: return prev + value;
        case OP_SUB: return prev - value;
        case OP_MAX: return max(value, prev);
        case OP_MIN: return min(value, prev);
        default:     return value;
    }
}

bool gto_op_reads_prev(uint operation) {
    return (operation & OP_ACCUM_MASK) != OP_SET;
}

// Squared distance from point p to the AABB [mn, mx].
float gto_aabb_dist2(float3 p, float3 mn, float3 mx) {
    float3 d = clamp(p, mn, mx) - p;
    return dot(d, d);
}

// ---------------------------------------------------------------------------
// Workgroup-wide exclusive scan of a uint2, wave-size agnostic.
// Requires GTO_SCAN_WG_SIZE to be defined to the workgroup size before inclusion.
// All threads of the group must call it (it contains barriers).
// ---------------------------------------------------------------------------
#ifdef GTO_SCAN_WG_SIZE
#define GTO_SCAN_MAX_WAVES 128
groupshared uint2 gto_scan_wave[GTO_SCAN_MAX_WAVES];
groupshared uint2 gto_scan_total;

uint2 gto_wg_exclusive_scan(uint2 v, uint tid, out uint2 total) {
    uint lane_count = WaveGetLaneCount();
    uint wave_idx   = tid / lane_count;
    uint2 wave_excl = uint2(WavePrefixSum(v.x), WavePrefixSum(v.y));
    uint2 wave_sum  = uint2(WaveActiveSum(v.x), WaveActiveSum(v.y));
    if (WaveIsFirstLane()) {
        gto_scan_wave[wave_idx] = wave_sum;
    }
    GroupMemoryBarrierWithGroupSync();
    if (tid == 0u) {
        uint num_waves = (GTO_SCAN_WG_SIZE + lane_count - 1u) / lane_count;
        uint2 run = uint2(0u, 0u);
        for (uint w = 0u; w < num_waves; ++w) {
            uint2 t = gto_scan_wave[w];
            gto_scan_wave[w] = run;
            run += t;
        }
        gto_scan_total = run;
    }
    GroupMemoryBarrierWithGroupSync();
    uint2 result = gto_scan_wave[wave_idx] + wave_excl;
    total = gto_scan_total;
    GroupMemoryBarrierWithGroupSync(); // gto_scan_wave may be reused by the next call
    return result;
}
#endif

#endif
