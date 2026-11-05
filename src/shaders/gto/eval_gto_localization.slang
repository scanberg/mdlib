// eval_gto_localization.slang
// The electron localization function (ELF) or the localized orbital locator (LOL) on a 3D grid, from the spins'
// density matrices over a contracted GTO basis (md_gto_gpu_localization_launch, md_gto.h).
//
// Per voxel, per spin s, from D_s and the AOs chi with their gradients:
//
//     rho_s      = sum_ij D_ij chi_i chi_j
//     grad rho_s = sum_ij D_ij (chi_i grad chi_j + grad chi_i chi_j)
//     tau_s      = 1/2 sum_ij D_ij grad chi_i . grad chi_j        (the positive definite kinetic energy density)
//
// then ELF or LOL as md_qm_elf / md_qm_lol take them, ELF faded to 0 where float runs out (md_qm_elf_noise) (the Kohout-Savin spin form; one matrix, the total density,
// for a closed shell, taken as two equal spins). The structure is the reference density kernel's
// (eval_gto_density.slang): the AOs that reach a group of 8^3 voxels screened, then tile pairs of 16 AOs, each
// AO's value and gradient evaluated in registers and contracted with the matrices' tiles in groupshared memory.

#include "../md_gpu.slang"
#include "gto_dmat_common.slang"

struct PGTO {
    float coeff;
    float alpha;
    uint ijkl;
};

struct RootArgs {
    float4x4 world_to_model;
    float4x4 index_to_world;
    float4   step;
    uint4    grid_dim;
    uint     num_cgtos;
    uint     num_spins;     // 1: D_a is the total density of a closed shell; 2: D_a and D_b are the spins'
    uint     lol;           // 0: ELF, 1: LOL
    float    rho_min;       // below it (both spins) the value is 0
    uint     dm_tiled;      // the matrices are tile pools of md_gto_gpu_dmats (dm_tao, dm_map of each)
    uint     dm_nt;
    float    fade_lo;       // ELF faded to 0 as its noise (localization_noise) goes from fade_lo to fade_hi; both 0: none
    float    fade_hi;

    uint*   cgto_atom_idx;
    float*  cgto_r;
    uint2*  cgto_off_len;
    PGTO*   pgto;
    float4* atom_xyz;
    float*  D_a;
    float*  D_b;
    uint*   dm_tao_a;
    uint*   dm_map_a;
    uint*   dm_tao_b;
    uint*   dm_map_b;
    DescriptorHandle<RWTexture3D<float>> out_tex;
};

struct MdGpuRoot { RootArgs* args_ptr; };
[[vk::push_constant]]
ConstantBuffer<MdGpuRoot> md_root;

#define args (*md_root.args_ptr)
#define out_vol (args.out_tex)

uint4 unpack_ijkl(uint packed) {
    return uint4((packed >> 0) & 0xFFu, (packed >> 8) & 0xFFu, (packed >> 16) & 0xFFu, (packed >> 24) & 0xFFu);
}

#define WG_X 8
#define WG_Y 8
#define WG_Z 8
#define WG_SIZE (WG_X * WG_Y * WG_Z)
#define TILE_SIZE 16
#define TILE_AREA (TILE_SIZE * TILE_SIZE)
#define MAX_SCREENED_CGTOS 2048
#define MAX_PGTOS_PER_CGTO 16
#define MAX_TILE_PGTOS (TILE_SIZE * MAX_PGTOS_PER_CGTO)

// 3/10 (6 pi^2)^(2/3): the uniform electron gas' kinetic energy density of one spin, over rho_s^(5/3)
#define ELF_CF_SPIN 4.557799872345596

groupshared uint  screened_cgtos[MAX_SCREENED_CGTOS];
groupshared uint  num_screened_cgtos;

groupshared float D_tile_a[TILE_SIZE][TILE_SIZE];
groupshared float D_tile_b[TILE_SIZE][TILE_SIZE];
groupshared float3 tile_cgto_xyz[TILE_SIZE];
groupshared uint2  tile_pgto_range[TILE_SIZE];

groupshared float pgtos_tile_coeff[MAX_TILE_PGTOS];
groupshared float pgtos_tile_alpha[MAX_TILE_PGTOS];
groupshared uint  pgtos_tile_ijkl [MAX_TILE_PGTOS];
groupshared uint  num_pgtos;
groupshared uint  pgtos_tile_offset[TILE_SIZE];

uint get_index(uint i, uint j, uint N) {
    uint row = (i < j) ? i : j;
    uint col = (i < j) ? j : i;
    return row * (2 * N - row + 1) / 2 + (col - row);
}

float d_at(float* D, uint* tao, uint* map, uint gi, uint gj) {
    return args.dm_tiled != 0u ? gto_dmat_at(tao, map, D, args.dm_nt, gi, gj) : D[get_index(gi, gj, args.num_cgtos)];
}

void fill_D_tiles(uint tile_i, uint tile_j, uint num_cgtos, uint tid) {
    uint baseI = tile_i * TILE_SIZE;
    uint baseJ = tile_j * TILE_SIZE;
    for (uint index = tid; index < TILE_AREA; index += WG_SIZE) {
        uint row = index / TILE_SIZE;
        uint col = index % TILE_SIZE;
        uint slot_i = baseI + row;
        uint slot_j = baseJ + col;
        float a = 0.0, b = 0.0;
        if (slot_i < num_cgtos && slot_j < num_cgtos) {
            // clamped: inactive lanes never form out-of-range addresses (llvmpipe evaluates their loads)
            uint gi = min(screened_cgtos[slot_i], args.num_cgtos - 1u);
            uint gj = min(screened_cgtos[slot_j], args.num_cgtos - 1u);
            a = d_at(args.D_a, args.dm_tao_a, args.dm_map_a, gi, gj);
            if (args.num_spins > 1u) b = d_at(args.D_b, args.dm_tao_b, args.dm_map_b, gi, gj);
        }
        D_tile_a[row][col] = a;
        D_tile_b[row][col] = b;
    }
}

// As the density kernel's: CGTO metadata and slots by the first TILE_SIZE threads, then the PGTOs by all
void fill_cgtos_tile(uint tile_offset, uint num_cgtos, uint tid) {
    if (tid < TILE_SIZE) {
        uint slot = tile_offset * TILE_SIZE + tid;
        uint src = 0u, len = 0u, dst = 0u;
        float3 center = float3(0.0f);
        if (slot < num_cgtos) {
            uint gi = min(screened_cgtos[slot], args.num_cgtos - 1u);
            uint2 off_len = args.cgto_off_len[gi];
            src = off_len.x;
            len = off_len.y;
            center = args.atom_xyz[args.cgto_atom_idx[gi]].xyz;
            InterlockedAdd(num_pgtos, len, dst);
        }
        pgtos_tile_offset[tid] = src;
        tile_pgto_range[tid] = uint2(dst, dst + len);
        tile_cgto_xyz[tid] = center;
    }
    GroupMemoryBarrierWithGroupSync();
    uint i = tid / MAX_PGTOS_PER_CGTO;
    if (i < TILE_SIZE) {
        uint j = tid % MAX_PGTOS_PER_CGTO;
        uint2 r = tile_pgto_range[i];
        if (r.x + j < r.y) {
            uint dst = r.x + j;
            PGTO p = args.pgto[pgtos_tile_offset[i] + j];
            pgtos_tile_coeff[dst] = p.coeff;
            pgtos_tile_alpha[dst] = p.alpha;
            pgtos_tile_ijkl [dst] = p.ijkl;
        }
    }
}

// x^n for n <= 5 (a monomial of a g function's derivative: one more)
float pow_le5(float x, uint n) {
    float x2 = x * x;
    switch (n) {
        case 0u: return 1.0f;
        case 1u: return x;
        case 2u: return x2;
        case 3u: return x2 * x;
        case 4u: return x2 * x2;
        default: return x2 * x2 * x;
    }
}

// d/dx of x^i e^(-a r^2), over e^(-a r^2): i x^(i-1) - 2 a x^(i+1)
float dmono(float x, uint i, float alpha) {
    float lo = i > 0u ? float(i) * pow_le5(x, i - 1u) : 0.0f;
    return lo - 2.0f * alpha * pow_le5(x, i + 1u);
}

// The value and the gradient of every AO of the tile at coord (0 for the slots past the screened ones)
void eval_phis(out float phi[TILE_SIZE], out float3 grad[TILE_SIZE], float3 coord) {
    for (uint i = 0u; i < TILE_SIZE; ++i) {
        float3 d  = coord - tile_cgto_xyz[i];
        float  r2 = dot(d, d);
        uint2 range = tile_pgto_range[i];
        uint  cnt   = range.y - range.x;
        float  v = 0.0f;
        float3 g = float3(0.0f);
        [unroll]
        for (uint j = 0u; j < MAX_PGTOS_PER_CGTO; ++j) {
            if (j >= cnt) break;
            uint  slot  = range.x + j;
            float alpha = pgtos_tile_alpha[slot];
            float e     = pgtos_tile_coeff[slot] * exp(-alpha * r2);
            uint4 ijk   = unpack_ijkl(pgtos_tile_ijkl[slot]);
            float fx = pow_le5(d.x, ijk.x), fy = pow_le5(d.y, ijk.y), fz = pow_le5(d.z, ijk.z);
            v += e * fx * fy * fz;
            g += e * float3(dmono(d.x, ijk.x, alpha) * fy * fz, fx * dmono(d.y, ijk.y, alpha) * fz, fx * fy * dmono(d.z, ijk.z, alpha));
        }
        phi[i]  = v;
        grad[i] = g;
    }
}

// md_qm_elf_noise and md_qm_elf_fade (md_qm.h): chi's error to expect, float's epsilon times tau / D0, and the
// weight that fades ELF out as it grows, 1 to 0 by a smoothstep in its log from lo to hi
#define ELF_NOISE_EPS 5.9604645e-8
float elf_fade(float noise, float lo, float hi) {
    if (!(noise < hi)) return 0.0f;     // NaN and inf too (rho^(5/3) underflowed)
    if (noise <= lo) return 1.0f;
    float t = (log(noise) - log(lo)) / (log(hi) - log(lo));
    return 1.0f - t * t * (3.0f - 2.0f * t);
}

// The terms of one spin, accumulated
struct Terms {
    float  rho;
    float3 grad;
    float  tau;
};

// Pair (i, j) once, of weight w (1 on the diagonal, 2 off it, as both (i, j) and (j, i)): of the ordered sums
// rho = sum D chi_i chi_j, grad rho = sum D (grad chi_i chi_j + chi_i grad chi_j), tau = 1/2 sum D grad chi_i . grad chi_j
void add_pair(inout Terms t, float D, float pi, float3 gi, float pj, float3 gj, float w) {
    float wd = w * D;
    t.rho  += wd * pi * pj;
    t.grad += wd * (pi * gj + pj * gi);
    t.tau  += (0.5f * wd) * dot(gi, gj);
}

[shader("compute")]
[numthreads(WG_X, WG_Y, WG_Z)]
void main(uint3 gid : SV_GroupID, uint3 gtid : SV_GroupThreadID, uint tid : SV_GroupIndex, uint3 dtid : SV_DispatchThreadID) {
    if (tid == 0) num_screened_cgtos = 0;
    GroupMemoryBarrierWithGroupSync();

    bool active = all(dtid < args.grid_dim.xyz);
    float3 model_aabb_min = float3( gid * uint3(WG_X, WG_Y, WG_Z))                 * args.step.xyz;
    float3 model_aabb_max = float3((gid + uint3(1,1,1)) * uint3(WG_X, WG_Y, WG_Z)) * args.step.xyz;

    // The AOs that reach the group's box (their radius: the basis's cutoff)
    for (uint i = tid; i < args.num_cgtos; i += WG_SIZE) {
        bool passed = false;
        if (i < args.num_cgtos) {
            float r = args.cgto_r[min(i, args.num_cgtos - 1u)];
            if (r != 0.0) {
                float3 center    = args.atom_xyz[args.cgto_atom_idx[min(i, args.num_cgtos - 1u)]].xyz;
                float3 model_xyz = mul(args.world_to_model, float4(center, 1.0)).xyz;
                float3 d         = clamp(model_xyz, model_aabb_min, model_aabb_max) - model_xyz;
                passed = dot(d, d) < r * r;
            }
        }
        uint wave_count  = WaveActiveCountBits(passed);
        uint wave_offset = WavePrefixSum(passed ? 1u : 0u);
        uint wave_base   = 0u;
        if (WaveIsFirstLane() && wave_count > 0u) InterlockedAdd(num_screened_cgtos, wave_count, wave_base);
        wave_base = WaveReadLaneFirst(wave_base);
        if (passed) {
            uint idx = wave_base + wave_offset;
            if (idx < MAX_SCREENED_CGTOS) screened_cgtos[idx] = i;
        }
    }
    GroupMemoryBarrierWithGroupSync();

    Terms ta = { 0.0f, float3(0.0f), 0.0f };
    Terms tb = { 0.0f, float3(0.0f), 0.0f };
    const bool two = args.num_spins > 1u;

    uint num_cgtos = min(num_screened_cgtos, MAX_SCREENED_CGTOS);
    if (num_cgtos > 0) {
        float  phi_mu[TILE_SIZE], phi_nu[TILE_SIZE];
        float3 grad_mu[TILE_SIZE], grad_nu[TILE_SIZE];
        float3 coord = mul(args.index_to_world, float4(float3(dtid), 1.0)).xyz;
        uint num_tiles = (num_cgtos + TILE_SIZE - 1) / TILE_SIZE;

        for (uint tile_i = 0; tile_i < num_tiles; ++tile_i) {
            // the diagonal tile: the pairs i <= j within it
            if (tid == 0) num_pgtos = 0;
            GroupMemoryBarrierWithGroupSync();
            fill_D_tiles(tile_i, tile_i, num_cgtos, tid);
            fill_cgtos_tile(tile_i, num_cgtos, tid);
            GroupMemoryBarrierWithGroupSync();
            if (active) {
                eval_phis(phi_mu, grad_mu, coord);
                for (uint i = 0; i < TILE_SIZE; ++i) {
                    add_pair(ta, D_tile_a[i][i], phi_mu[i], grad_mu[i], phi_mu[i], grad_mu[i], 1.0f);
                    if (two) add_pair(tb, D_tile_b[i][i], phi_mu[i], grad_mu[i], phi_mu[i], grad_mu[i], 1.0f);
                    for (uint j = i + 1; j < TILE_SIZE; ++j) {
                        add_pair(ta, D_tile_a[i][j], phi_mu[i], grad_mu[i], phi_mu[j], grad_mu[j], 2.0f);
                        if (two) add_pair(tb, D_tile_b[i][j], phi_mu[i], grad_mu[i], phi_mu[j], grad_mu[j], 2.0f);
                    }
                }
            }
            // the tiles after it: every pair
            for (uint tile_j = tile_i + 1; tile_j < num_tiles; ++tile_j) {
                GroupMemoryBarrierWithGroupSync();   // the tiles are read until here
                if (tid == 0) num_pgtos = 0;
                GroupMemoryBarrierWithGroupSync();
                fill_D_tiles(tile_i, tile_j, num_cgtos, tid);
                fill_cgtos_tile(tile_j, num_cgtos, tid);
                GroupMemoryBarrierWithGroupSync();
                if (active) {
                    eval_phis(phi_nu, grad_nu, coord);
                    for (uint i = 0; i < TILE_SIZE; ++i) {
                        for (uint j = 0; j < TILE_SIZE; ++j) {
                            add_pair(ta, D_tile_a[i][j], phi_mu[i], grad_mu[i], phi_nu[j], grad_nu[j], 2.0f);
                            if (two) add_pair(tb, D_tile_b[i][j], phi_mu[i], grad_mu[i], phi_nu[j], grad_nu[j], 2.0f);
                        }
                    }
                }
            }
            GroupMemoryBarrierWithGroupSync();   // before the next diagonal tile's fill
        }
    }

    if (active) {
        // A closed shell: two equal spins of half of it (rho / 2, |grad rho|^2 / 4, tau / 2)
        float ra, g2a, tau_a, rb, g2b, tau_b;
        if (two) {
            ra = ta.rho; g2a = dot(ta.grad, ta.grad); tau_a = ta.tau;
            rb = tb.rho; g2b = dot(tb.grad, tb.grad); tau_b = tb.tau;
        } else {
            ra = 0.5f * ta.rho; g2a = 0.25f * dot(ta.grad, ta.grad); tau_a = 0.5f * ta.tau;
            rb = ra; g2b = g2a; tau_b = tau_a;
        }
        float value = 0.0f;
        float rho = ra + rb;
        if (rho > 0.0f && rho >= args.rho_min) {
            float d0 = float(ELF_CF_SPIN) * ((ra > 0.0f ? pow(ra, 5.0f / 3.0f) : 0.0f) + (rb > 0.0f ? pow(rb, 5.0f / 3.0f) : 0.0f));
            if (args.lol != 0u) {
                float tau = tau_a + tau_b;
                float t = tau > 0.0f ? d0 / tau : 0.0f;
                value = tau > 0.0f ? t / (1.0f + t) : 1.0f;
            } else if (d0 > 0.0f) {
                float da = ra > 0.0f ? max(tau_a - g2a / (8.0f * ra), 0.0f) : 0.0f;
                float db = rb > 0.0f ? max(tau_b - g2b / (8.0f * rb), 0.0f) : 0.0f;
                float chi = (da + db) / d0;
                value = 1.0f / (1.0f + chi * chi);
                if (args.fade_hi > 0.0f) value *= elf_fade(ELF_NOISE_EPS * (tau_a + tau_b) / d0, args.fade_lo, args.fade_hi);
            }
        }
        out_vol[dtid] = value;
    }
}
