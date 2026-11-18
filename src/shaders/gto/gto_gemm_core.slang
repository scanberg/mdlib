// gto_gemm_core.slang
//
// THE GEMM CORE of the two-pass grid evaluations (count + Phi passes, eval_gto_density_gemm.slang): one main loop,
// specialised at compile time by what it multiplies and what is made of the product.
//
//   X[r][m] = sum_k A[r][k] B[k][m]     r: the GP rows of a group, m < ncols the columns, k < nk (the block's AOs)
//
// The product goes through the columns in tiles of GM, and for each tile through k in chunks of GK held in groupshared
// memory, the next chunk loaded into registers while the current one is multiplied. TPX x TPY threads (16 x 16 unless
// the kernel says otherwise); thread (tx, ty) holds the RP x RM part of X of rows ty * RP + i and columns
// m0 + tx * RM + j, and reads its rows and columns of the chunk as float4 where RP, RM are multiples of 4. After each column tile the EPILOGUE
// takes that part of X; it may read anything (A at the same rows, say) and keeps what it accumulates in its own state.
//
//   OPERANDS   IGemmOperands: A[r][k] and B[k][m] as the loads read them (consecutive threads read consecutive rows of
//              A and consecutive columns of B: both streams coalesced where the operands are laid out so). The core
//              keeps the loads in range (k < nk, m < ncols) and, for a symmetric product, to k >= m only.
//   EPILOGUE   IGemmEpilogue: the part of X of each column tile.
//
// What uses it: the density (eval_gto_density_core.slang: symmetric, ncols = nk), ELF and LOL
// (eval_gto_localization_core.slang: symmetric, four rows a voxel), the orbitals (eval_gto_mo_gemm.slang: B the
// orbitals' coefficients, ncols the orbitals), a density matrix's tiles from orbitals (gto_dmat.slang: P_I^T diag(w) P_J,
// 16 x 128 per group of 64 threads).
//
// A SYMMETRIC product (B = w o D of a symmetric D, w 1 on the diagonal and 2 off it, k >= m only) sums
// sum_{m <= k} w D[k][m] A[r][k] - half the multiply-adds - and is exact for a quadratic form a^T D a = sum_m a[m] X[m]
// with X of the same a; for a symmetric bilinear form 2 a^T D b = sum_m (a[m] Xb[m] + b[m] Xa[m]) with the X of both.
//
// What it costs: GP x GM x n / 2 (symmetric) multiply-adds per column tile, of which the epilogue's are a few per
// element of X. The specialisations (Slang generics, instantiated per kernel) are each a kernel of their own, with the
// operands' and the epilogue's code inlined into the loop: one source, as many binaries as before.
//
// Before including: GP (rows per group, a multiple of TPY), GM (columns per tile, a multiple of TPX), optionally TPX and
// TPY (16 each by default). The kernel is numthreads(TPX, TPY, 1), and passes its SV_GroupIndex and SV_GroupThreadID.xy.
// groupshared memory is reached through gc_ld / gc_st / gc_ld4 (gc_smem* directly is for the core alone).

#ifndef GTO_GEMM_CORE_SLANG
#define GTO_GEMM_CORE_SLANG

// For measurements: bumped with every change to the core's loop
#define GTO_GEMM_CORE_REVISION 3

#define GK   16                 // k per chunk
#ifndef TPX
#define TPX  16                 // threads along the columns
#endif
#ifndef TPY
#define TPY  16                 // threads along the rows
#endif
#define NT   (TPX * TPY)
#define RP   (GP / TPY)         // rows per thread
#define RM   (GM / TPX)         // columns per thread
#define PF_A ((GK * GP) / NT)   // A chunk elements per thread
#define PF_B ((GK * GM) / NT)   // B chunk elements per thread

// [0, GK*GP)            : the A chunk, [kk][r]
// [GK*GP, GK*(GP+GM))   : the B chunk, [kk][mm]
// after the core        : the epilogue's reductions (GK * (GP + GM) floats)
// Written a float at a time, read as float4 where a thread's rows (columns) allow. SPIR-V: float4 elements; Metal:
// floats, since Metal may write one component of a threadgroup vector as a read-modify-write of the whole vector, and
// two threads writing components of one float4 would then lose a write (as gto_iso_gpu.slang's K_BUF).
#define GC_FLOATS (GK * GP + GK * GM)
groupshared float4 gc_smem4[GC_FLOATS / 4];
groupshared float  gc_smem1[GC_FLOATS];
float gc_ld(uint i) {
    __target_switch { case spirv: return gc_smem4[i >> 2][i & 3u];
                      default:    return gc_smem1[i]; }
}
float4 gc_ld4(uint i4) {
    __target_switch { case spirv: return gc_smem4[i4];
                      default:    return float4(gc_smem1[4u * i4], gc_smem1[4u * i4 + 1u], gc_smem1[4u * i4 + 2u], gc_smem1[4u * i4 + 3u]); }
}
void gc_st(uint i, float v) {
    __target_switch { case spirv: gc_smem4[i >> 2][i & 3u] = v;
                      default:    gc_smem1[i] = v; }
}

// A thread's rows (columns) of a chunk read as float4: RP (RM) a multiple of 4 (the chunk's offsets are then too)
#if (RP % 4) == 0
#define GC_A_VEC 1
#else
#define GC_A_VEC 0
#endif
#if (RM % 4) == 0
#define GC_B_VEC 1
#else
#define GC_B_VEC 0
#endif

interface IGemmOperands {
    // A[row][k], row < GP, k < nk
    float a(uint row, uint k);
    // B[k][m], k < nk, m < ncols (symmetric: k >= m, weighted)
    float b(uint k, uint m);
};

interface IGemmEpilogue {
    // The part of X of the columns [m0, m0 + GM): X[i][j] at row ty * RP + i, column m0 + tx * RM + j (only those
    // below ncols are of the product)
    [mutating] void take(uint m0, uint ncols, uint tx, uint ty, float X[RP][RM]);
};

[ForceInline]
void gemm_core<O : IGemmOperands, E : IGemmEpilogue>(O ops, inout E epi, uint nk, uint ncols, bool symmetric, uint tid, uint tx, uint ty) {
    for (uint m0 = 0u; m0 < ncols; m0 += GM) {
        float X[RP][RM];
        [ForceUnroll] for (int i = 0; i < RP; ++i)
            [ForceUnroll] for (int j = 0; j < RM; ++j) X[i][j] = 0.0f;

        // the next chunk, in registers
        float pf_a[PF_A];
        float pf_b[PF_B];

        // Loads the chunk of k from k0 into pf_a / pf_b: consecutive threads, consecutive rows (A) and columns (B)
        #define GC_LOAD_CHUNK(k0_)                                                                  \
        {                                                                                           \
            [ForceUnroll] for (int r = 0; r < PF_A; ++r) {                                          \
                uint idx = tid + uint(r) * NT;                                                      \
                uint k   = (k0_) + idx / GP;                                                        \
                float v  = ops.a(idx % GP, min(k, nk - 1u));                                        \
                pf_a[r]  = (k < nk) ? v : 0.0f;                                                     \
            }                                                                                       \
            [ForceUnroll] for (int r = 0; r < PF_B; ++r) {                                          \
                uint idx = tid + uint(r) * NT;                                                      \
                uint k   = (k0_) + idx / GM;                                                        \
                uint m   = m0 + idx % GM;                                                           \
                float v  = 0.0f;                                                                    \
                if (k < nk && m < ncols && (!symmetric || k >= m)) v = ops.b(k, m);                 \
                pf_b[r] = v;                                                                        \
            }                                                                                       \
        }

        // a symmetric product starts at the tile's own rows of B: the ones above are zero
        const uint k_first = symmetric ? m0 : 0u;
        if (k_first < nk) {                      // nk 0 (an empty block of the orbitals): nothing to load
            GC_LOAD_CHUNK(k_first);
        }

        for (uint k0 = k_first; k0 < nk; k0 += GK) {
            GroupMemoryBarrierWithGroupSync();   // the previous chunk consumed
            [ForceUnroll] for (int r = 0; r < PF_A; ++r) gc_st(tid + uint(r) * NT, pf_a[r]);
            [ForceUnroll] for (int r = 0; r < PF_B; ++r) gc_st(GK * GP + tid + uint(r) * NT, pf_b[r]);
            GroupMemoryBarrierWithGroupSync();

            if (k0 + GK < nk) {
                GC_LOAD_CHUNK(k0 + GK);           // in flight while this chunk is multiplied
            }

            [ForceUnroll]
            for (int kk = 0; kk < GK; ++kk) {
                float a[RP];
                float bb[RM];
#if GC_A_VEC
                [ForceUnroll] for (int q = 0; q < RP / 4; ++q) {
                    const float4 v = gc_ld4((uint(kk) * GP + ty * RP) / 4u + uint(q));
                    a[4 * q] = v.x; a[4 * q + 1] = v.y; a[4 * q + 2] = v.z; a[4 * q + 3] = v.w;
                }
#else
                [ForceUnroll] for (int i = 0; i < RP; ++i) a[i]  = gc_ld(uint(kk) * GP + ty * RP + uint(i));
#endif
#if GC_B_VEC
                [ForceUnroll] for (int q = 0; q < RM / 4; ++q) {
                    const float4 v = gc_ld4((GK * GP + uint(kk) * GM + tx * RM) / 4u + uint(q));
                    bb[4 * q] = v.x; bb[4 * q + 1] = v.y; bb[4 * q + 2] = v.z; bb[4 * q + 3] = v.w;
                }
#else
                [ForceUnroll] for (int j = 0; j < RM; ++j) bb[j] = gc_ld(GK * GP + uint(kk) * GM + tx * RM + uint(j));
#endif
                [ForceUnroll] for (int i = 0; i < RP; ++i)
                    [ForceUnroll] for (int j = 0; j < RM; ++j) X[i][j] += a[i] * bb[j];
            }
        }
        #undef GC_LOAD_CHUNK

        epi.take(m0, ncols, tx, ty, X);
    }
}

// The sum over the TPX threads of each row of one value per row (v[i]: row ty * RP + i), into out_row (row r's sum,
// read by thread r < GP after this returns; every thread calls it). Through gc_smem, so after the core.
[ForceInline]
void gemm_row_sums(float v[RP], uint tid, uint tx, uint ty, out float out_row) {
    GroupMemoryBarrierWithGroupSync();   // gc_smem free
    [ForceUnroll] for (int i = 0; i < RP; ++i) gc_st((ty * RP + uint(i)) * TPX + tx, v[i]);
    GroupMemoryBarrierWithGroupSync();
    out_row = 0.0f;
    if (tid < GP) {
        [ForceUnroll] for (int k = 0; k < TPX; ++k) out_row += gc_ld(tid * TPX + uint(k));
    }
}

#endif
