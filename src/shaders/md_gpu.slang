// md_gpu.slang
//
// Prelude for compute shaders launched through md_gpu.h.
//
// There is deliberately almost nothing here. Textures and samplers are passed
// as Slang `DescriptorHandle<T>` values living directly in the argument
// struct, so there is no binding table to declare, nothing to keep in sync
// with the host, and no combinatorial list of (dimension x sampled type x
// access) arrays. Slang emits whatever descriptor arrays a given shader needs
// and md_gpu provides a single heap for them to index.
//
// Usage:
//
//     #include "md_gpu.slang"
//
//     struct Args {
//         uint4  dim;                  // xyz used, w spare (see below)
//         DescriptorHandle<RWTexture3D<float>> vol;
//         DescriptorHandle<Texture3D<float4>>  src;
//         DescriptorHandle<SamplerState>       samp;
//         float* dst;
//     };
//     MD_KERNEL_ARGS(Args);
//
//     [shader("compute")]
//     [numthreads(4, 4, 4)]
//     void main(uint3 tid : SV_DispatchThreadID) {
//         Args   a   = MD_ARGS;
//         uint3  dim = a.dim.xyz;
//         if (any(tid >= dim)) return;
//         a.vol[tid] = a.src.SampleLevel(a.samp, float3(tid) / float3(dim), 0).x;
//     }
//
// A handle is type-checked: a DescriptorHandle<RWTexture3D<float>> cannot be
// used where a Texture2D is expected, so the mismatched-index class of bug
// that an untyped `uint` handle allows is a compile error here.
//
// Build requirement: shaders must be compiled with
//
//     -bindless-space-index MD_GPU_BINDLESS_SPACE   (currently 0)
//
// which is what cmake/CompileGpuShaders.cmake passes. Combined with the
// getDescriptorFromHandle override below, Slang then places its sampler heap
// at (space, binding 0), sampled images at binding 2 and storage images at
// binding 3; md_gpu creates a matching descriptor set. See md_gpu_vulkan.c
// "Bindless heap" for the host side of that contract.

#ifndef MD_GPU_SLANG_INCLUDED
#define MD_GPU_SLANG_INCLUDED

// ---------------------------------------------------------------------------
// Argument block
// ---------------------------------------------------------------------------
//
// A kernel takes exactly one pointer to a caller-defined struct. md_gpu
// copies the struct into device memory and passes its address in an 8-byte
// push constant -- the only push constant used, and small enough to fit the
// 128-byte Vulkan minimum with room to spare.
//
// MIRRORING THE ARGUMENT STRUCT IN C
//
// The two backends do not lay vectors out identically. Slang emits SPIR-V with
// *scalar* layout -- a vector's alignment is its component's -- while MSL gives
// vec2 8/8 and vec4 16/16. The sizes agree; only the alignment differs. So one
// C struct is correct on both backends exactly when every vector member sits at
// an offset that satisfies the stricter (Metal) alignment:
//
//     vec2, texture/sampler handle    offset must be a multiple of 8
//     vec4, float4x4                  offset must be a multiple of 16
//
// Mirror those members with md_gpu.h's md_gpu_uint4 / md_gpu_float4x4 / ...
// which declare that alignment, and the C struct falls out correct with no
// per-backend #ifdef anywhere:
//
//     struct Args {                        // Slang
//         uint4  dim;                      // xyz used
//         float  scale;
//         uint   _pad;
//         DescriptorHandle<RWTexture3D<float>> vol;
//         float* dst;
//     };
//
//     typedef struct {                     // C
//         md_gpu_uint4 dim;
//         float        scale;
//         uint32_t     _pad;
//         md_gpu_tex_t vol;
//         uint64_t     dst;
//     } my_args_t;
//
// Scalars and device pointers place themselves. Texture and sampler handles do
// not: Slang lowers DescriptorHandle<T> to two 32-bit words, which SPIR-V
// aligns to 4 and MSL to 8, so a handle after an odd number of 32-bit scalars
// lands 4 bytes apart on the two targets. Pad it to an 8-byte offset, exactly
// as for a vec2.
//
// NO THREE-COMPONENT VECTORS IN ARGUMENT STRUCTS. uint3 is 12 bytes in SPIR-V
// and 16 in MSL, and padding cannot close that: MSL's slack sits *inside* the
// vector and SPIR-V's sits after it, so every following member shifts by 4 on
// one target only -- putting the vector on a 16-byte boundary fixes where it
// starts and nothing else:
//
//     struct { uint*a; uint*b; uint3 v; uint pad; uint next; }
//         SPIR-V: v@16  pad@28  next@32
//         MSL:    v@16  pad@32  next@36
//
// Use uint4 and ignore .w, or three scalars. Everywhere else -- locals,
// groupshared, arithmetic -- 3-vectors are fine and cost nothing.
//
// This is enforced, not trusted: tools/check_gpu_arg_layout.py compiles each
// kernel for both targets and rejects any argument struct containing a
// 3-vector, or whose two layouts disagree at any member. If it passes, one C
// mirror written with the md_gpu_* types is correct on every backend.

// ---------------------------------------------------------------------------
// Bindless binding assignment
// ---------------------------------------------------------------------------
//
// Slang's default SPIR-V lowering for DescriptorHandle<T> is
// BindlessDescriptorOptions.VkMutable: every resource type is aliased onto a
// single binding, which is only legal with VK_EXT_mutable_descriptor_type.
// That extension is unavailable on NVIDIA hardware before Turing, and Khronos
// explicitly advises against mutable descriptors when the alternative is a
// handful of concrete types rather than emulating a DirectX 12 heap -- static
// types let implementations keep their fast paths.
//
// So take the other preset. `None` gives one binding per descriptor type:
//
//     0  SAMPLER            3  STORAGE_IMAGE          6  UNIFORM_BUFFER
//     1  COMBINED_SAMPLER   4  UNIFORM_TEXEL_BUFFER   7  STORAGE_BUFFER
//     2  SAMPLED_IMAGE      5  STORAGE_TEXEL_BUFFER   8  other
//
// Nothing else changes: handles stay typed, the argument structs are
// untouched, and Metal is unaffected (this option is SPIR-V only, and the
// override lowers to an identity function there).
export T getDescriptorFromHandle<T>(DescriptorHandle<T> handle) where T : IOpaqueDescriptor
{
    return defaultGetDescriptorFromHandle(handle, BindlessDescriptorOptions.None);
}

// Declares the push constant holding the pointer to this kernel's argument
// struct. Place once at file scope, after defining the struct.
#define MD_KERNEL_ARGS(T)                                    \
    struct MdGpuRoot { T* args; };                           \
    [[vk::push_constant]] ConstantBuffer<MdGpuRoot> md_root

// The argument struct, by value. Load it once into a local.
#define MD_ARGS (*md_root.args)

// The argument struct, by pointer, when you want to write through it.
#define MD_ARGS_PTR (md_root.args)

#endif // MD_GPU_SLANG_INCLUDED
