// eval_gto_mo_shell.slang
//
// Molecular orbitals on a 3D grid, one workgroup per 8x8x8 block (32x16x1 for a grid of one layer):
//
//   psi_m(r) = sum_mu C[m][mu] phi_mu(r)
//   result   = sum_m w_m * psi_m          (mode 0, PSI)
//            = sum_m w_m * psi_m^2        (mode 1, PSI_SQUARED)
//
// with w_m = 1 unless a weight buffer is given (weights may be negative, e.g. for
// differences of orbital densities; they are applied after squaring).
//
// Compared with the per-AO reference kernel it replaced (eval_gto_mo.slang, gone):
//   * Shell based: screening per contracted shell, the radial part evaluated once per
//     shell and shared by its Cartesian components.
//   * No per-tile barriers. Shells are screened in chunks of one workgroup; the
//     passing shells of a chunk go to groupshared memory and are evaluated straight
//     away, reading their primitives and coefficients with workgroup-uniform loads.
//     Two barriers per chunk of up to WG_SIZE shells instead of four per 16 AOs, and
//     no limit on the number of shells per block (no silent truncation).
//   * Coefficient-aware screening (optional): a shell is skipped where
//     max|C| * |phi| < cutoff rather than |phi| < cutoff, using the shell's most
//     diffuse exponent to shrink (or grow) its radius:  r'^2 = r^2 + ln(max|C|) / alpha_min.
//     Localised orbitals then touch far fewer shells per block.
//   * VPT voxels per thread share every shell setup and coefficient load.
//   * MB orbitals per pass (1 or 8): MB = 1 for the common single orbital without
//     wasted accumulators; MB = 8 for small sets of orbitals. Larger sets belong on
//     the GEMM path (eval_gto_mo_gemm.slang).
//
// GRAD (the *_grad.slang wrappers): the gradient of the result too, into out_grad (xyz, w 0), by the same op -
//   grad psi_m = sum_mu C[m][mu] grad phi_mu, grad phi = N (grad(x^i y^j z^k) R + x^i y^j z^k 2 (r - A) dR/dr^2),
//   dR/dr^2 = -sum_p c_p alpha_p exp(-alpha_p r^2); grad (psi^2) = 2 psi grad psi. World units (per bohr).

#include "../md_gpu.slang"

#ifndef VPT
#define VPT 2
#endif
#ifndef MB
#define MB 1
#endif

// A group's block of voxels: 8x8x8, VPT voxels a thread along z; or for a grid of one layer (the *_2d.slang
// wrappers) 32x16x1, VPT along y. A kernel per shape: one for both, the shape an argument as the GEMM path's
// (gto_gemm_args.slang), was 5% behind for one orbital on an M4 Pro's volume and 9-19% for 8 on a GTX 1060.
#ifndef BLOCK_X
#define BLOCK_X 8
#define BLOCK_Y 8
#define BLOCK_Z 8
#endif
#define BLOCK_SHAPE uint3(BLOCK_X, BLOCK_Y, BLOCK_Z)
#if BLOCK_Z > 1
#define WG_X BLOCK_X
#define WG_Y BLOCK_Y
#define WG_Z (BLOCK_Z / VPT)
#define VPT_STEP uint3(0u, 0u, uint(WG_Z))
#else
#define WG_X BLOCK_X
#define WG_Y (BLOCK_Y / VPT)
#define WG_Z 1
#define VPT_STEP uint3(0u, uint(WG_Y), 0u)
#endif
#define WG_SIZE (WG_X * WG_Y * WG_Z)

#define GTO_SCAN_WG_SIZE WG_SIZE
#include "gto_common.slang"

struct RootArgs {
    float4x4 world_to_model;
    float4x4 index_to_world;
    float4   step;
    uint4    grid_dim;
    uint     num_shells;
    uint     num_aos;
    uint     num_mos;
    uint     mode;              // 0 = psi, 1 = psi^2
    uint     operation;
    uint     has_weights;
    uint     coeff_screening;
    uint     _pad0;

    uint*    shell_atom_idx;
    float*   shell_r;
    uint4*   shell_info;        // (prim_offset, num_prims, l, ao_offset)
    float2*  shell_prim;        // (coeff, alpha * log2(e))
    float4*  atom_xyz;
    float*   coeffs;            // [num_mos][num_aos]
    float*   weights;           // [num_mos], if has_weights
    DescriptorHandle<RWTexture3D<float>> out_tex;
    DescriptorHandle<RWTexture3D<float4>> out_grad;     // GRAD only
};

struct MdGpuRoot { RootArgs* args_ptr; };
[[vk::push_constant]]
ConstantBuffer<MdGpuRoot> md_root;

#define args (*md_root.args_ptr)

// Cartesian components in the order of md_gto.h, packed i | j << 4 | k << 8, and their
// normalisation 1/sqrt((2i-1)!!(2j-1)!!(2k-1)!!). Component c of a shell with angular
// momentum l is entry CART_BASE[l] + c.
static const uint CART_BASE[5] = { 0u, 1u, 4u, 10u, 20u };
static const uint CART_IJK[35] = {
    0x000u, 0x001u, 0x010u, 0x100u, 0x002u, 0x011u, 0x101u, 0x020u, 0x110u, 0x200u, 0x003u, 0x012u,
    0x102u, 0x021u, 0x111u, 0x201u, 0x030u, 0x120u, 0x210u, 0x300u, 0x004u, 0x013u, 0x103u, 0x022u,
    0x112u, 0x202u, 0x031u, 0x121u, 0x211u, 0x301u, 0x040u, 0x130u, 0x220u, 0x310u, 0x400u };
static const float CART_NORM[35] = {
    1.000000000f, 1.000000000f, 1.000000000f, 1.000000000f, 0.577350269f, 1.000000000f, 1.000000000f,
    0.577350269f, 1.000000000f, 0.577350269f, 0.258198890f, 0.577350269f, 0.577350269f, 0.577350269f,
    1.000000000f, 0.577350269f, 0.258198890f, 0.577350269f, 0.577350269f, 0.258198890f, 0.097590007f,
    0.258198890f, 0.258198890f, 0.333333333f, 0.577350269f, 0.333333333f, 0.258198890f, 0.577350269f,
    0.577350269f, 0.258198890f, 0.097590007f, 0.258198890f, 0.333333333f, 0.258198890f, 0.097590007f };

#define LOG2E 1.44269504089f
#define LN2   0.69314718056f

#ifdef GRAD
// The gradient of x^i y^j z^k (ijk packed as CART_IJK)
float3 gto_monomial_grad(float3 d, uint ijk) {
    const uint i = ijk & 0xFu, j = (ijk >> 4u) & 0xFu, k = (ijk >> 8u) & 0xFu;
    return float3(i != 0u ? float(i) * gto_monomial(d, ijk - 0x001u) : 0.0f,
                  j != 0u ? float(j) * gto_monomial(d, ijk - 0x010u) : 0.0f,
                  k != 0u ? float(k) * gto_monomial(d, ijk - 0x100u) : 0.0f);
}
#endif

groupshared uint chunk_shell[WG_SIZE];

[shader("compute")]
[numthreads(WG_X, WG_Y, WG_Z)]
void main(uint3 gid : SV_GroupID, uint3 gtid : SV_GroupThreadID, uint tid : SV_GroupIndex) {
    const float3 aabb_min = float3( gid             * BLOCK_SHAPE) * args.step.xyz;
    const float3 aabb_max = float3((gid + uint3(1u)) * BLOCK_SHAPE) * args.step.xyz;

    float3 coord[VPT];
    uint3  vox[VPT];
    float  result[VPT];
    [ForceUnroll] for (int v = 0; v < VPT; ++v) {
        vox[v]    = gid * BLOCK_SHAPE + gtid + uint(v) * VPT_STEP;
        coord[v]  = mul(args.index_to_world, float4(float3(vox[v]), 1.0f)).xyz;
        result[v] = 0.0f;
    }
#ifdef GRAD
    float3 result_g[VPT];
    [ForceUnroll] for (int v = 0; v < VPT; ++v) result_g[v] = float3(0.0f);
#endif

    for (uint mo_base = 0u; mo_base < args.num_mos; mo_base += MB) {
        const uint mo_count = min((uint)MB, args.num_mos - mo_base);

        float psi[VPT][MB];
        [ForceUnroll] for (int v = 0; v < VPT; ++v)
            [ForceUnroll] for (int m = 0; m < MB; ++m) psi[v][m] = 0.0f;
#ifdef GRAD
        float3 gpsi[VPT][MB];
        [ForceUnroll] for (int v = 0; v < VPT; ++v)
            [ForceUnroll] for (int m = 0; m < MB; ++m) gpsi[v][m] = float3(0.0f);
#endif

        for (uint base = 0u; base < args.num_shells; base += WG_SIZE) {
            // Screen one chunk of shells.
            // Loads use a clamped index so that no lane, active or not, ever forms an
            // out-of-range address (some implementations evaluate loads of inactive lanes).
            const uint s_raw = base + tid;
            const uint s = min(s_raw, args.num_shells - 1u);
            bool pass = false;
            if (s_raw < args.num_shells) {
                const float r = args.shell_r[s];
                if (r > 0.0f) {
                    float r2 = r * r;
                    if (args.coeff_screening != 0u) {
                        const uint4 info = args.shell_info[s];
                        const uint  nc   = gto_num_cart(info.z);
                        float cmax = 0.0f;
                        for (uint m = 0u; m < mo_count; ++m) {
                            const uint row = (mo_base + m) * args.num_aos + info.w;
                            for (uint c = 0u; c < nc; ++c) cmax = max(cmax, abs(args.coeffs[row + c]));
                        }
                        float amin = 3.4e38f;
                        for (uint p = 0u; p < info.y; ++p) amin = min(amin, args.shell_prim[info.x + p].y);
                        // log(0) = -inf removes shells whose coefficients are all zero.
                        r2 += log(cmax) * (LOG2E / amin);
                    }
                    const float3 c = args.atom_xyz[args.shell_atom_idx[s]].xyz;
                    const float3 m = mul(args.world_to_model, float4(c, 1.0f)).xyz;
                    pass = gto_aabb_dist2(m, aabb_min, aabb_max) < r2;
                }
            }
            uint2 total;
            const uint2 excl = gto_wg_exclusive_scan(uint2(pass ? 1u : 0u, 0u), tid, total);
            if (pass) chunk_shell[excl.x] = s;   // s == s_raw for passing lanes
            GroupMemoryBarrierWithGroupSync();

            // Evaluate the chunk's shells. Every load below is uniform across the group.
            for (uint k = 0u; k < total.x; ++k) {
                // Clamped like s above: always in range, so a speculatively evaluated
                // dependent load cannot fault (seen with llvmpipe).
                const uint   sh   = min(WaveReadLaneFirst(chunk_shell[k]), args.num_shells - 1u);
                const uint4  info = args.shell_info[sh];
                const float3 ctr  = args.atom_xyz[args.shell_atom_idx[sh]].xyz;

                float3 d[VPT];
                float  r2[VPT];
                float  R[VPT];
                [ForceUnroll] for (int v = 0; v < VPT; ++v) {
                    d[v]  = coord[v] - ctr;
                    r2[v] = dot(d[v], d[v]);
                    R[v]  = 0.0f;
                }
#ifdef GRAD
                float Rp[VPT];      // sum_p c_p alpha_p log2(e) exp(-alpha_p r^2): dR/dr^2 = -Rp ln 2
                [ForceUnroll] for (int v = 0; v < VPT; ++v) Rp[v] = 0.0f;
                for (uint p = 0u; p < info.y; ++p) {
                    const float2 cp = args.shell_prim[info.x + p];
                    [ForceUnroll] for (int v = 0; v < VPT; ++v) {
                        const float e = cp.x * exp2(-cp.y * r2[v]);
                        R[v]  += e;
                        Rp[v] += cp.y * e;
                    }
                }
#else
                for (uint p = 0u; p < info.y; ++p) {
                    const float2 cp = args.shell_prim[info.x + p];
                    [ForceUnroll] for (int v = 0; v < VPT; ++v) R[v] += cp.x * exp2(-cp.y * r2[v]);
                }
#endif

                const uint nc = gto_num_cart(info.z);
                const uint tb = CART_BASE[min(info.z, 4u)];
                for (uint c = 0u; c < nc; ++c) {
                    const uint  ijk = CART_IJK[tb + c];
                    const float nrm = CART_NORM[tb + c];
                    float cf[MB];
                    [ForceUnroll] for (int m = 0; m < MB; ++m) {
                        // In-range row for every m (the select below may evaluate both sides).
                        const float cv = args.coeffs[(mo_base + min(uint(m), mo_count - 1u)) * args.num_aos + info.w + c];
                        cf[m] = (uint(m) < mo_count) ? cv : 0.0f;
                    }
                    [ForceUnroll] for (int v = 0; v < VPT; ++v) {
                        const float mono = gto_monomial(d[v], ijk);
                        const float phi = nrm * mono * R[v];
                        [ForceUnroll] for (int m = 0; m < MB; ++m) psi[v][m] += cf[m] * phi;
#ifdef GRAD
                        const float3 gphi = nrm * (gto_monomial_grad(d[v], ijk) * R[v] - (2.0f * LN2 * mono * Rp[v]) * d[v]);
                        [ForceUnroll] for (int m = 0; m < MB; ++m) gpsi[v][m] += cf[m] * gphi;
#endif
                    }
                }
            }
            GroupMemoryBarrierWithGroupSync();   // chunk_shell is rewritten by the next chunk
        }

        [ForceUnroll] for (int m = 0; m < MB; ++m) {
            if (uint(m) < mo_count) {
                float w = 1.0f;
                if (args.has_weights != 0u) w = args.weights[mo_base + uint(m)];
                [ForceUnroll] for (int v = 0; v < VPT; ++v) {
                    const float x = psi[v][m];
                    result[v] += w * ((args.mode != 0u) ? x * x : x);
#ifdef GRAD
                    result_g[v] += w * ((args.mode != 0u) ? 2.0f * x * gpsi[v][m] : gpsi[v][m]);
#endif
                }
            }
        }
    }

    [ForceUnroll] for (int v = 0; v < VPT; ++v) {
        if (all(vox[v] < args.grid_dim.xyz)) {
            const float prev = gto_op_reads_prev(args.operation) ? args.out_tex[vox[v]] : 0.0f;
            args.out_tex[vox[v]] = gto_apply_op(result[v], prev, args.operation);
#ifdef GRAD
            const float3 prev_g = gto_op_reads_prev(args.operation) ? args.out_grad[vox[v]].xyz : float3(0.0f);
            args.out_grad[vox[v]] = float4(gto_apply_op_grad(result[v], result_g[v], prev, prev_g, args.operation), 0.0f);
#endif
        }
    }
}
