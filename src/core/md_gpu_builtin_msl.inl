/* Metal Shading Language source for the built-in md_gpu_make_grid helper.
   Compiled at device creation with newLibraryWithSource:. The Vulkan backend
   uses the SPIR-V in md_gpu_builtin_spv.inl, generated from
   src/shaders/gpu/md_gpu_make_grid.slang; this is the hand-written MSL
   equivalent of that same kernel. */

static const char* md_gpu_make_grid_msl =
"#include <metal_stdlib>\n"
"using namespace metal;\n"
"\n"
"struct MdMakeGridArgs {\n"
"    device uint* count;\n"
"    device uint* out_grid;\n"
"    uint3        local;\n"
"    uint         _pad;\n"
"};\n"
"\n"
"kernel void md_gpu_make_grid(constant MdMakeGridArgs& args [[buffer(0)]]) {\n"
"    uint n  = args.count[0];\n"
"    uint lx = max(args.local.x, 1u);\n"
"    uint ly = max(args.local.y, 1u);\n"
"    uint lz = max(args.local.z, 1u);\n"
"    args.out_grid[0] = (n + lx - 1u) / lx;\n"
"    args.out_grid[1] = (1u + ly - 1u) / ly;\n"
"    args.out_grid[2] = (1u + lz - 1u) / lz;\n"
"}\n";
