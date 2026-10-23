// gto_dmat.slang
// A density matrix in 16 x 16 tiles of the AO pairs that meet (md_gto_gpu_dmat, md_gto.c), built from
// orbitals: tile (I, J) = P_I^T diag(w) P_J.
//
// P holds the orbitals in the tiled AO order, P[k * ntao + t] for orbital k and tiled AO t (ntao = 16 nt,
// zero past the last AO). A work item is a row tile I and up to 8 of its column tiles (consecutive
// entries of the row's list): a 16 x 128 block. A group of 64 threads does one, each thread a 4 x 8 block
// (rows 4 rq .. 4 rq + 3, columns 8 co .. 8 co + 7 of the 128), so per orbital a thread makes 32
// multiply-adds from 3 groupshared loads. The orbitals go through groupshared memory KC at a time, the
// next chunk loaded into registers while the current one is used. Float accumulation.
//
// r2. (r1: 256 threads, 2 x 4 per thread, 2 groupshared loads per 8 multiply-adds: bound by groupshared
// bandwidth, 108 G multiply-adds/s on a GTX 1060.)

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

#define KC 16
groupshared float4 sa[KC * 4];          // [k][4]: the row tile's 16 AOs, weighted
groupshared float4 sb[KC * 32];         // [k][32]: the 8 column tiles' 128 AOs
groupshared uint   s_cols[8];           // the column tile of each slot

// orbital k's 4 coefficients of tiled AOs t .. t + 3
float4 dm_load4(float* P, uint k, uint ntao, uint t) {
    float* src = P + k * ntao + t;
    return float4(src[0], src[1], src[2], src[3]);
}

[shader("compute")]
[numthreads(64, 1, 1)]
void dmat_tiles(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) {
    DmatArgs a = *dmat_root.args;
    uint item = gid.y * a.items_x + gid.x;
    if (item >= a.num_items) return;                       // uniform per group
    uint4 it = a.items[item];
    uint I = it.x, e0 = it.y, cnt = it.z;
    uint r = a.num_mo;
    if (tid < 8u) s_cols[tid] = tid < cnt ? a.cols[e0 + tid] : 0u;
    GroupMemoryBarrierWithGroupSync();

    uint co = tid & 15u, rq = tid >> 4;                    // columns 8 co .., rows 4 rq ..
    uint slot = co >> 1;

    // the loads of a chunk: sa by (k = tid / 4, quad tid % 4); sb by i = tid + 64 j: (k = i / 32, slot (i / 4) % 8, quad i % 4)
    float4 pa;
    float4 pb[8];
    #define DM_LOAD(k0_)                                                                        \
    {                                                                                           \
        uint kn_ = min(uint(KC), r - (k0_));                                                    \
        uint ka_ = tid >> 2;                                                                    \
        uint kk_ = min((k0_) + ka_, r - 1u);    /* clamped: no lane forms an address past P */  \
        float4 va_ = dm_load4(a.P, kk_, a.ntao, I * 16u + (tid & 3u) * 4u) * a.w[kk_];          \
        pa = ka_ < kn_ ? va_ : float4(0.0);                                                     \
        [ForceUnroll] for (int j_ = 0; j_ < 8; ++j_) {                                          \
            uint i_ = tid + 64u * uint(j_);                                                     \
            uint kb_ = i_ >> 5, s_ = (i_ >> 2) & 7u;                                            \
            float4 vb_ = dm_load4(a.P, min((k0_) + kb_, r - 1u), a.ntao, s_cols[s_] * 16u + (i_ & 3u) * 4u); \
            pb[j_] = (kb_ < kn_ && s_ < cnt) ? vb_ : float4(0.0);                               \
        }                                                                                       \
    }

    float4 acc[4][2];
    [ForceUnroll] for (int i = 0; i < 4; ++i) { acc[i][0] = float4(0.0); acc[i][1] = float4(0.0); }

    DM_LOAD(0u);
    for (uint k0 = 0u; k0 < r; k0 += uint(KC)) {
        GroupMemoryBarrierWithGroupSync();                 // the previous chunk used
        sa[tid] = pa;
        [ForceUnroll] for (int j = 0; j < 8; ++j) sb[tid + 64u * uint(j)] = pb[j];
        GroupMemoryBarrierWithGroupSync();
        if (k0 + uint(KC) < r) DM_LOAD(k0 + uint(KC));      // in flight while this chunk is used
        uint kn = min(uint(KC), r - k0);
        for (uint k = 0u; k < kn; ++k) {
            float4 x = sa[k * 4u + rq];
            float4 b0 = sb[k * 32u + co * 2u];
            float4 b1 = sb[k * 32u + co * 2u + 1u];
            [ForceUnroll] for (int i = 0; i < 4; ++i) {
                acc[i][0] += x[i] * b0;
                acc[i][1] += x[i] * b1;
            }
        }
    }
    #undef DM_LOAD

    if (slot < cnt) {
        float* dst = a.tiles + (e0 + slot) * 256u + (rq * 4u) * 16u + (co & 1u) * 8u;
        [ForceUnroll] for (int i = 0; i < 4; ++i) {
            float* row = dst + uint(i) * 16u;
            row[0] = acc[i][0].x; row[1] = acc[i][0].y; row[2] = acc[i][0].z; row[3] = acc[i][0].w;
            row[4] = acc[i][1].x; row[5] = acc[i][1].y; row[6] = acc[i][1].z; row[7] = acc[i][1].w;
        }
    }
}
