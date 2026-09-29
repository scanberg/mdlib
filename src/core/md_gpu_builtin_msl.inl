/* Metal Shading Language source for the built-in md_gpu_make_grid helper and
   the byte-granular copy/fill kernel (md_gpu_byte_op, Metal only).
   Compiled at device creation with newLibraryWithSource:. The Vulkan backend
   uses the SPIR-V in md_gpu_builtin_spv.inl, generated from
   src/shaders/gpu/md_gpu_make_grid.slang; this is the hand-written MSL
   equivalent of that same kernel.

   Being hand-written, it has to imitate Slang's lowering rather than inherit
   it: the bound buffer holds a *pointer* to the argument struct, matching what
   MD_KERNEL_ARGS produces for every other kernel. It previously took the struct
   by reference, which quietly made this the one kernel the old
   struct-at-buffer-0 binding was correct for.

   `local` is a uint4 with an unused .w rather than a uint3, because argument
   structs may not contain 3-vectors -- SPIR-V would give one 12 bytes and MSL
   16, and nothing after it would line up. Keep this in step with Args in
   src/shaders/gpu/md_gpu_make_grid.slang by hand; only the SPIR-V half is
   generated. */

static const char* md_gpu_make_grid_msl =
"#include <metal_stdlib>\n"
"using namespace metal;\n"
"\n"
"struct MdMakeGridArgs {\n"
"    device uint* count;\n"
"    device uint* out_grid;\n"
"    uint4        local;\n"
"};\n"
"\n"
"struct MdGpuRoot { device MdMakeGridArgs* args; };\n"
"\n"
"kernel void md_gpu_make_grid(constant MdGpuRoot& root [[buffer(0)]]) {\n"
"    device MdMakeGridArgs* args = root.args;\n"
"    uint n  = args->count[0];\n"
"    uint lx = max(args->local.x, 1u);\n"
"    uint ly = max(args->local.y, 1u);\n"
"    uint lz = max(args->local.z, 1u);\n"
"    args->out_grid[0] = (n + lx - 1u) / lx;\n"
"    args->out_grid[1] = (1u + ly - 1u) / ly;\n"
"    args->out_grid[2] = (1u + lz - 1u) / lz;\n"
"}\n";

/* Byte copy / fill for what blits may not do on macOS: offsets or sizes that
   are not multiples of 4. Each thread handles 16 bytes (MD_MTL_BYTE_OP_SPAN);
   src == 0 means fill with `value`. Keep in step with md_mtl_byte_op_args_t. */
static const char* md_gpu_byte_op_msl =
"#include <metal_stdlib>\n"
"using namespace metal;\n"
"\n"
"struct MdByteOpArgs {\n"
"    device uchar*       dst;\n"
"    device const uchar* src;\n"
"    uint                size;\n"
"    uint                value;\n"
"};\n"
"\n"
"struct MdGpuRoot { device MdByteOpArgs* args; };\n"
"\n"
"kernel void md_gpu_byte_op(constant MdGpuRoot& root [[buffer(0)]],\n"
"                           uint tid [[thread_position_in_grid]]) {\n"
"    device MdByteOpArgs* a = root.args;\n"
"    uint begin = tid * 16u;\n"
"    if (begin >= a->size) return;\n"
"    uint end = min(begin + 16u, a->size);\n"
"    if (a->src) {\n"
"        for (uint i = begin; i < end; ++i) a->dst[i] = a->src[i];\n"
"    } else {\n"
"        uchar v = uchar(a->value);\n"
"        for (uint i = begin; i < end; ++i) a->dst[i] = v;\n"
"    }\n"
"}\n";
