// gto_dmat.slang
// A density matrix in 16 x 16 tiles of the AO pairs that meet (md_gto_gpu_dmat, md_gto.c), built from
// orbitals: tile (I, J) = P_I^T diag(w) P_J.
//
// P holds the orbitals in the tiled AO order, P[k * ntao + t] for orbital k and tiled AO t (ntao = 16 nt,
// zero past the last AO). A work item is a row tile I and up to 8 of its column tiles (consecutive
// entries of the row's list): a 16 x 128 block, on the GEMM core (gto_gemm_core.slang): a group of 16 x 4
// threads, each a 4 x 8 block, the orbitals through groupshared memory 16 at a time, the next chunk loaded
// into registers while the current one is used. Float accumulation.
//
// r3. (r1: 256 threads, 2 x 4 per thread: bound by groupshared bandwidth, 108 G multiply-adds/s on a GTX 1060.
// r2: its own loop, 64 threads of 4 x 8 with float4 groupshared reads; the core's r3 the same bits, 9.71 ms against
// 9.21 on a protein on an RTX 3080, and r2 retired.)

#include "../md_gpu.slang"

struct DmatArgs {
    uint   num_mo;
    uint   num_items;
    uint   items_x;      // groups along x: item = gid.y * items_x + gid.x
    uint   ntao;
    float* P;
    float* w;
    uint4* items;        // x: row tile, y: first entry, z: entries (1 .. 8)
    uint*  cols;         // column tile of each entry
    float* tiles;        // entry e: tiles[e * 256 + 16 a + c]
};

struct DmatRoot { DmatArgs* args; };
[[vk::push_constant]] ConstantBuffer<DmatRoot> dmat_root;

// ---------------------------------------------------------------------------------------------------------------------
// The work items on the GEMM core. Rows: the row tile's 16 AOs (A = w o P_I^T), columns: the item's up to 8
// column tiles, 16 AOs each (B = P_J), k: the orbitals.
// ---------------------------------------------------------------------------------------------------------------------
#define GP  16
#define GM  128
#include "gto_gemm_core.slang"

#define GC_TX 16
#define GC_TY 4
#define RP (GP / GC_TY)
#define RM (GM / GC_TX)
GC_MEM(GcMem, GP, GM)

#define dargs (*dmat_root.args)
groupshared uint dc_cols[8];        // the column tile of each of the item's slots

struct DmatOperands : IGemmOperands {
    uint I;
    float a(uint row, uint k) {
        return dargs.w[k] * dargs.P[k * dargs.ntao + I * 16u + row];
    }
    float b(uint k, uint m) {
        return dargs.P[k * dargs.ntao + dc_cols[m >> 4] * 16u + (m & 15u)];
    }
};

// X is the item's 16 x (16 cnt) block: tile (I, J of slot) entry e0 + slot, [a][c] at 16 a + c
struct DmatEpilogue : IGemmEpilogue {
    uint e0;
    [mutating] void take<let R : int, let M : int>(uint m0, uint ncols, uint tx, uint ty, float X[R][M]) {
        [ForceUnroll] for (int j = 0; j < M; ++j) {
            const uint m = m0 + tx * M + uint(j);
            if (m < ncols) {
                float* dst = dargs.tiles + (e0 + (m >> 4)) * 256u + (m & 15u);
                [ForceUnroll] for (int i = 0; i < R; ++i) dst[(ty * R + uint(i)) * 16u] = X[i][j];
            }
        }
    }
};

[shader("compute")]
[numthreads(GC_TX, GC_TY, 1)]
void dmat_tiles_core(uint3 gid : SV_GroupID, uint3 gtid : SV_GroupThreadID, uint tid : SV_GroupIndex) {
    const uint item = gid.y * dargs.items_x + gid.x;
    if (item >= dargs.num_items) return;                    // uniform per group
    const uint4 it = dargs.items[item];
    const uint cnt = it.z;
    if (tid < 8u) dc_cols[tid] = tid < cnt ? dargs.cols[it.y + tid] : 0u;
    GroupMemoryBarrierWithGroupSync();

    DmatOperands ops;
    ops.I = it.x;
    DmatEpilogue epi;
    epi.e0 = it.y;
    gemm_core<GP, GM, GC_TX, GC_TY, GcMem, DmatOperands, DmatEpilogue>(ops, epi, dargs.num_mo, 16u * cnt, false, tid, gtid.x, gtid.y);
}
