// gto_gemm_core.slang
//
// THE GEMM CORE of the two-pass grid evaluations (count + Phi passes, eval_gto_density_gemm.slang): one main loop,
// specialised at compile time by what it multiplies and what is made of the product.
//
//   X[r][m] = sum_k A[r][k] B[k][m]     r: the GP rows of a group, m < ncols the columns, k < nk (the block's AOs)
//
// The product goes through the columns in tiles of GM, and for each tile through k in chunks of GK held in groupshared
// memory, the next chunk loaded into registers while the current one is multiplied. 16 x 16 threads; thread (tx, ty)
// holds the RP x RM part of X of rows ty * RP + i and columns m0 + tx * RM + j. After each column tile the EPILOGUE
// takes that part of X; it may read anything (A at the same rows, say) and keeps what it accumulates in its own state.
//
//   OPERANDS   IGemmOperands: A[r][k] and B[k][m] as the loads read them (consecutive threads read consecutive rows of
//              A and consecutive columns of B: both streams coalesced where the operands are laid out so). The core
//              keeps the loads in range (k < nk, m < ncols) and, for a symmetric product, to k >= m only.
//   EPILOGUE   IGemmEpilogue: the part of X of each column tile.
//
// What uses it: the density (eval_gto_density_core.slang: symmetric, ncols = nk), ELF and LOL
// (eval_gto_localization_core.slang: symmetric, four rows a voxel), the orbitals (eval_gto_mo_gemm.slang: B the
// orbitals' coefficients, ncols the orbitals).
//
// A SYMMETRIC product (B = w o D of a symmetric D, w 1 on the diagonal and 2 off it, k >= m only) sums
// sum_{m <= k} w D[k][m] A[r][k] - half the multiply-adds - and is exact for a quadratic form a^T D a = sum_m a[m] X[m]
// with X of the same a; for a symmetric bilinear form 2 a^T D b = sum_m (a[m] Xb[m] + b[m] Xa[m]) with the X of both.
//
// What it costs: GP x GM x n / 2 (symmetric) multiply-adds per column tile, of which the epilogue's are a few per
// element of X. The specialisations (Slang generics, instantiated per kernel) are each a kernel of their own, with the
// operands' and the epilogue's code inlined into the loop: one source, as many binaries as before.
//
// Before including: GP (rows per group, a multiple of 16 dividing what the kernel's rows are), GM (columns per tile, a
// multiple of 16), and gto_gemm_args.slang. The kernel is numthreads(16, 16, 1).

#ifndef GTO_GEMM_CORE_SLANG
#define GTO_GEMM_CORE_SLANG

// For measurements: bumped with every change to the core's loop
#define GTO_GEMM_CORE_REVISION 2

#define GK   16                 // k per chunk
#define TPX  16                 // threads along the columns
#define TPY  16                 // threads along the rows
#define NT   (TPX * TPY)
#define RP   (GP / TPY)         // rows per thread
#define RM   (GM / TPX)         // columns per thread
#define PF_A ((GK * GP) / NT)   // A chunk elements per thread
#define PF_B ((GK * GM) / NT)   // B chunk elements per thread

// [0, GK*GP)            : the A chunk, [kk][r]
// [GK*GP, GK*(GP+GM))   : the B chunk, [kk][mm]
// after the core        : the epilogue's reductions (GK * (GP + GM) floats)
groupshared float gc_smem[GK * GP + GK * GM];

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
            [ForceUnroll] for (int r = 0; r < PF_A; ++r) gc_smem[tid + uint(r) * NT] = pf_a[r];
            [ForceUnroll] for (int r = 0; r < PF_B; ++r) gc_smem[GK * GP + tid + uint(r) * NT] = pf_b[r];
            GroupMemoryBarrierWithGroupSync();

            if (k0 + GK < nk) {
                GC_LOAD_CHUNK(k0 + GK);           // in flight while this chunk is multiplied
            }

            [ForceUnroll]
            for (int kk = 0; kk < GK; ++kk) {
                float a[RP];
                float bb[RM];
                [ForceUnroll] for (int i = 0; i < RP; ++i) a[i]  = gc_smem[uint(kk) * GP + ty * RP + uint(i)];
                [ForceUnroll] for (int j = 0; j < RM; ++j) bb[j] = gc_smem[GK * GP + uint(kk) * GM + tx * RM + uint(j)];
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
    [ForceUnroll] for (int i = 0; i < RP; ++i) gc_smem[(ty * RP + uint(i)) * TPX + tx] = v[i];
    GroupMemoryBarrierWithGroupSync();
    out_row = 0.0f;
    if (tid < GP) {
        [ForceUnroll] for (int k = 0; k < TPX; ++k) out_row += gc_smem[tid * TPX + uint(k)];
    }
}

#endif
