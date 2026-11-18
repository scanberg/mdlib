// The block shape an argument (BLOCK_RUNTIME, gto_gemm_args.slang): one kernel for 8x8x8 blocks and a grid of one
// layer's 32x16x1. See eval_gto_density_gemm.slang.
#define BLOCK_RUNTIME 1
#include "eval_gto_density_gemm.slang"
