// md_gpu_make_grid.slang
//
// Built-in helper backing md_gpu_make_grid(). Reads a single uint32 element
// count from device memory and writes an indirect dispatch grid:
//
//     out[0] = ceil(count / local.x)
//     out[1] = ceil(1     / local.y)      (== 1 for the usual case)
//     out[2] = ceil(1     / local.z)
//
// This exists so that turning a device-side count into an indirect launch does
// not require every caller to write the same three-line shader.

#include "../md_gpu.slang"

struct Args {
    uint*  count;      // in : single uint32 element count
    uint*  out_grid;   // out: 3 x uint32 dispatch dimensions
    uint4  local;      // xyz = threads per group; w unused -- argument structs
                       // may not contain 3-vectors (see md_gpu.slang)
};
MD_KERNEL_ARGS(Args);

[shader("compute")]
[numthreads(1, 1, 1)]
void main() {
    Args a = MD_ARGS;
    uint n = a.count[0];
    uint lx = max(a.local.x, 1u);
    uint ly = max(a.local.y, 1u);
    uint lz = max(a.local.z, 1u);
    a.out_grid[0] = (n + lx - 1u) / lx;
    a.out_grid[1] = (1u + ly - 1u) / ly;
    a.out_grid[2] = (1u + lz - 1u) / lz;
}
