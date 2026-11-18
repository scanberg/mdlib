// eval_gto_localization_core.slang
//
// ELF and LOL on the GEMM core (gto_gemm_core.slang): per spin s, with D = D_s and the AOs phi with their gradients,
//
//   rho      = sum_mn D phi_m phi_n
//   grad rho = sum_mn D (grad phi_m phi_n + phi_m grad phi_n)
//   tau      = 1/2 sum_mn D grad phi_m . grad phi_n
//
// as one symmetric product (half the multiply-adds) over four rows per voxel: G_c = phi (c = 0) and d phi / d x_c
// (c = 1, 2, 3), all four from the Phi pass with gradients (phi_grad_main: four rows of Phi per AO),
//
//   X_c[p][m] = sum_{n >= m} w D[n][m] G_c[p][n]              (w: 1 on the diagonal, 2 off it)
//   rho       = sum_m G_0 X_0
//   d rho/dx_c = sum_m (G_c X_0 + G_0 X_c)                       (a symmetric bilinear form, from the X of both)
//   tau       = 1/2 sum_c sum_m G_c X_c
//
// 4 n^2 / 2 multiply-adds a voxel per spin: four times the density's, the least a gradient of tau's kind takes. The
// three terms of each voxel and spin go to the batch's terms (args.terms), and finish_main makes ELF or LOL of them
// (md_qm_elf, md_qm_lol: the Kohout-Savin spin form, as eval_gto_localization.slang).
//
// The rows of a group: GP of them, four channels of GP / 4 voxels, channel-major (row r: channel r / (GP / 4), voxel
// r % (GP / 4)), so that the loads of A read runs of GP / 4 voxels of one row of Phi. A thread's RP rows are of one
// channel (RP divides GP / 4).
//
// Compile-time configuration (the eval_gto_localization_core_*.slang wrappers): GP rows per group, GM columns per tile.

#include "gto_gemm_args.slang"

#ifndef GP
#define GP 128
#endif
#ifndef GM
#define GM 64
#endif

#include "gto_gemm_core.slang"

#define GC_TX 16
#define GC_TY 16
#define RP (GP / GC_TY)
#define RM (GM / GC_TX)
GC_MEM(GcMem, GP, GM)

#define LC_VOX (GP / 4)     // voxels per group
#define ELF_CF_SPIN 4.557799872345596   // 3/10 (6 pi^2)^(2/3): the uniform electron gas' kinetic energy density of one spin over rho_s^(5/3)

// Phi with gradients: row c of AO k at voxel p of the block
float lc_phi(uint phi_off, uint k, uint c, uint p) {
    return args.phi[phi_off + (k * 4u + c) * BLOCK_VOXELS + p];
}

// A: the four channels of the group's voxels; B: w o D_b of the spin
struct LocalizationOperands : IGemmOperands {
    uint phi_off, db_off, n, p0;
    float a(uint row, uint k) {
        return lc_phi(phi_off, k, row / LC_VOX, p0 + row % LC_VOX);
    }
    float b(uint k, uint m) {
        const float w = (k == m) ? 1.0f : 2.0f;
        return w * args.D_blocks[args.d_spin_offset + db_off + k * n + m];
    }
};

// Per row: channel 0 - rho and the X_0 halves of the three derivatives (sum G_c X_0); channel c > 0 - the X_c half of
// its derivative (sum G_0 X_c) and its part of 2 tau (sum G_c X_c)
struct LocalizationEpilogue : IGemmEpilogue {
    float acc[RP][4];
    uint  phi_off, p0;
    [mutating] void take<let R : int, let M : int>(uint m0, uint n, uint tx, uint ty, float X[R][M]) {
        const uint c = (ty * R) / LC_VOX;      // the thread's channel: uniform per wave of 32 or more
        [ForceUnroll] for (int j = 0; j < M; ++j) {
            uint m = m0 + tx * M + uint(j);
            if (m < n) {
                [ForceUnroll] for (int i = 0; i < R; ++i) {
                    const uint  p = p0 + (ty * R + uint(i)) % LC_VOX;
                    const float x = X[i][j];
                    if (c == 0u) {
                        acc[i][0] += lc_phi(phi_off, m, 0u, p) * x;
                        acc[i][1] += lc_phi(phi_off, m, 1u, p) * x;
                        acc[i][2] += lc_phi(phi_off, m, 2u, p) * x;
                        acc[i][3] += lc_phi(phi_off, m, 3u, p) * x;
                    } else {
                        acc[i][0] += lc_phi(phi_off, m, 0u, p) * x;
                        acc[i][1] += lc_phi(phi_off, m, c, p) * x;
                    }
                }
            }
        }
    }
};

[shader("compute")]
[numthreads(16, 16, 1)]
void main(uint3 gid : SV_GroupID, uint3 gtid : SV_GroupThreadID, uint tid : SV_GroupIndex) {
    const uint groups_per_block = BLOCK_VOXELS / LC_VOX;
    const uint b  = gid.x / groups_per_block;
    const uint pt = gid.x % groups_per_block;
    if (b >= args.num_blocks) return;   // uniform per group

    uint4 e0 = args.block_table[2u * b + 0u];
    uint4 e1 = args.block_table[2u * b + 1u];
    const uint n     = e0.w;
    const uint entry = e1.w;            // within the batch: where its terms go

    LocalizationOperands ops;
    ops.phi_off = e1.x;
    ops.db_off  = e1.z;
    ops.n       = n;
    ops.p0      = pt * LC_VOX;
    LocalizationEpilogue epi;
    epi.phi_off = e1.x;
    epi.p0      = pt * LC_VOX;
    [ForceUnroll] for (int i = 0; i < RP; ++i) [ForceUnroll] for (int q = 0; q < 4; ++q) epi.acc[i][q] = 0.0f;

    gemm_core<GP, GM, GC_TX, GC_TY, GcMem, LocalizationOperands, LocalizationEpilogue>(ops, epi, n, n, true, tid, gtid.x, gtid.y);

    // the rows' sums, four of each, then the voxels' terms of them
    float s[4];
    [ForceUnroll] for (int q = 0; q < 4; ++q) {
        float v[RP];
        [ForceUnroll] for (int i = 0; i < RP; ++i) v[i] = epi.acc[i][q];
        gemm_row_sums<GP, GC_TX, RP, GcMem>(v, tid, gtid.x, gtid.y, s[q]);
    }
    GroupMemoryBarrierWithGroupSync();
    const uint RS = GP * GC_TX;           // past what gemm_row_sums used
    if (tid < GP) [ForceUnroll] for (int q = 0; q < 4; ++q) GcMem::st(RS + tid * 4u + uint(q), s[q]);
    GroupMemoryBarrierWithGroupSync();
    if (tid < LC_VOX) {
        const uint v = tid;
        #define LC_RS(c, q) GcMem::ld(RS + ((c) * LC_VOX + v) * 4u + (q))
        const float rho = LC_RS(0u, 0u);
        const float gx  = LC_RS(0u, 1u) + LC_RS(1u, 0u);
        const float gy  = LC_RS(0u, 2u) + LC_RS(2u, 0u);
        const float gz  = LC_RS(0u, 3u) + LC_RS(3u, 0u);
        const float tau = 0.5f * (LC_RS(1u, 1u) + LC_RS(2u, 1u) + LC_RS(3u, 1u));
        #undef LC_RS
        const uint plane = args.terms_stride * BLOCK_VOXELS;
        const uint o = (args.spin * 5u) * plane + entry * BLOCK_VOXELS + pt * LC_VOX + v;
        args.terms[o]             = rho;
        args.terms[o + plane]     = gx;
        args.terms[o + 2u * plane] = gy;
        args.terms[o + 3u * plane] = gz;
        args.terms[o + 4u * plane] = tau;
    }
}

// ELF or LOL of the terms of a batch's entries (or 0 for the empty blocks: terms_valid 0), into out_tex: one group per
// block
[shader("compute")]
[numthreads(BLOCK_X, BLOCK_Y, BLOCK_Z)]
void finish_main(uint3 gid : SV_GroupID, uint3 gtid : SV_GroupThreadID, uint tid : SV_GroupIndex) {
    const uint b = gid.x;
    if (b >= args.num_blocks) return;
    uint4 e0 = args.block_table[2u * b + 0u];
    uint4 e1 = args.block_table[2u * b + 1u];
    const uint3 vox = e0.xyz * BLOCK_SHAPE + gtid;
    if (!all(vox < args.grid_dim.xyz)) return;
    float value = 0.0f;
    if (args.terms_valid != 0u) {
        const uint plane = args.terms_stride * BLOCK_VOXELS;
        const uint o = e1.w * BLOCK_VOXELS + tid;
        float ra = args.terms[o], tau_a = args.terms[o + 4u * plane];
        float3 ga = float3(args.terms[o + plane], args.terms[o + 2u * plane], args.terms[o + 3u * plane]);
        float g2a = dot(ga, ga);
        float rb, g2b, tau_b;
        if (args.num_spins > 1u) {
            const uint ob = o + 5u * plane;
            rb = args.terms[ob];
            tau_b = args.terms[ob + 4u * plane];
            float3 gb = float3(args.terms[ob + plane], args.terms[ob + 2u * plane], args.terms[ob + 3u * plane]);
            g2b = dot(gb, gb);
        } else {
            // a closed shell: two equal spins of half of it (rho / 2, |grad rho|^2 / 4, tau / 2)
            ra *= 0.5f; g2a *= 0.25f; tau_a *= 0.5f;
            rb = ra; g2b = g2a; tau_b = tau_a;
        }
        const float rho = ra + rb;
        if (rho > 0.0f && rho >= args.rho_min) {
            const float d0 = float(ELF_CF_SPIN) * ((ra > 0.0f ? pow(ra, 5.0f / 3.0f) : 0.0f) + (rb > 0.0f ? pow(rb, 5.0f / 3.0f) : 0.0f));
            if (args.lol != 0u) {
                const float tau = tau_a + tau_b;
                const float t = tau > 0.0f ? d0 / tau : 0.0f;
                value = tau > 0.0f ? t / (1.0f + t) : 1.0f;
            } else if (d0 > 0.0f) {
                const float da = ra > 0.0f ? max(tau_a - g2a / (8.0f * ra), 0.0f) : 0.0f;
                const float db = rb > 0.0f ? max(tau_b - g2b / (8.0f * rb), 0.0f) : 0.0f;
                const float chi = (da + db + args.d_min) / d0;
                value = 1.0f / (1.0f + chi * chi);     // (chi^2 inf where d0 underflows: 0)
            }
        }
    }
    args.out_tex[vox] = value;
}
