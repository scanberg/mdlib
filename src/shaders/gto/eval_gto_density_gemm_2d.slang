// Configuration wrapper, see eval_gto_density_gemm.slang: blocks of 32x16x1 voxels, for a grid of one layer.
#define BLOCK_X 32
#define BLOCK_Y 16
#define BLOCK_Z 1
#include "eval_gto_density_gemm.slang"
