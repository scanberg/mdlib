// md_gpu.slang
//
// Prelude for kernels and vertex/fragment shaders used through md_gpu.h.
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
//     MD_SHADER_ARGS(Args);
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

// Vectors through pointers: Slang's natural layout aligns a vector to its
// element, so a float4's (or uint4's) load or store through a pointer, alone
// or as a struct's member, is marked Aligned 4, and NVIDIA's compiler splits
// it into words (it may not assume more). loadAligned<16>(p) and
// storeAligned<16>(p, v) say more where it is known. Measure before using
// them: on a GTX 1060 the frames' gather in md_render_mol went from 6.3 to
// 2.8 ms with its store aligned, while md_render's view and bounds aligned
// gained nothing and took more registers (a 128-bit load wants four in a row,
// even for a kernel that uses three of them).

// Declares the push constant holding the pointer to the argument struct.
// Place once at file scope, after defining the struct. Every entry point in
// the file -- kernel, vertex or fragment -- reads the same root, which is how
// a draw's vertex and fragment stages share one argument struct.
#define MD_SHADER_ARGS(T)                                    \
    struct MdGpuRoot { T* args; };                           \
    [[vk::push_constant]] ConstantBuffer<MdGpuRoot> md_root

// The previous name, kept until callers have moved.
#define MD_KERNEL_ARGS(T) MD_SHADER_ARGS(T)

// The argument struct, by value. Load it once into a local.
#define MD_ARGS (*md_root.args)

// The argument struct, by pointer, when you want to write through it.
#define MD_ARGS_PTR (md_root.args)

// ---------------------------------------------------------------------------
// Draw-relative vertex and instance ids
// ---------------------------------------------------------------------------
//
// Direct draws (md_gpu_draw / md_gpu_draw_indexed) have no base vertex or base
// instance, so SV_VertexID and SV_InstanceID are already draw-relative on
// every backend. Indirect commands carry bases, and there the targets differ:
// Slang's SPIR-V subtracts them, Metal's [[vertex_id]] / [[instance_id]]
// include them. Shaders drawn indirectly pass the bases in and get the
// draw-relative id everywhere:
//
//     VOut vs_main(uint vid : SV_VertexID, uint iid : SV_InstanceID,
//                  uint base_vertex : SV_StartVertexLocation,
//                  uint base_instance : SV_StartInstanceLocation) {
//         uint v    = md_draw_vertex(vid, base_vertex);
//         uint i    = md_draw_instance(iid, base_instance);
//         Draw draw = MD_ARGS.draws[base_instance];      // one record per command
//         ...
//
// SV_StartInstanceLocation itself is portable, which makes it the way to give
// each indirect command its own record: a culling kernel writes command i with
// first_instance = i.
uint md_draw_vertex(uint vid, uint base_vertex) {
    __target_switch {
    case metal: return vid - base_vertex;
    default:    return vid;
    }
}

uint md_draw_instance(uint iid, uint base_instance) {
    __target_switch {
    case metal: return iid - base_instance;
    default:    return iid;
    }
}

// Half floats as bits in a uint, both ways, exactly (round to nearest even;
// subnormals, infinities and NaNs kept), in 32-bit arithmetic alone. Slang's
// f16tof32 and f32tof16 go through 16-bit types, which make a SPIR-V module
// declare the Float16 and Int16 capabilities: a device must enable
// shaderFloat16 and shaderInt16 for it, which not every one has (a GTX 1060
// reports neither enabled; Vulkan validation flags every such module).
float md_half_to_float(uint h) {
    uint s = (h & 0x8000u) << 16u;
    uint e = (h >> 10u) & 0x1Fu;
    uint m = h & 0x3FFu;
    if (e == 0u) {                              // zero and subnormals: m 2^-24
        float v = float(m) * 5.9604644775390625e-8;
        return s != 0u ? -v : v;
    }
    if (e == 31u) return asfloat(s | 0x7F800000u | (m << 13u));
    return asfloat(s | ((e + 112u) << 23u) | (m << 13u));
}

uint md_float_to_half(float f) {
    uint x = asuint(f);
    uint s = (x >> 16u) & 0x8000u;
    uint e = (x >> 23u) & 0xFFu;
    uint m = x & 0x7FFFFFu;
    if (e == 0xFFu) return s | 0x7C00u | (m != 0u ? 0x200u : 0u);
    int he = int(e) - 112;
    if (he >= 31) return s | 0x7C00u;           // too large: infinity
    if (he <= 0) {                              // a subnormal, or zero
        if (he < -10) return s;
        m |= 0x800000u;
        uint sh = uint(14 - he);
        uint hm = m >> sh, rem = m & ((1u << sh) - 1u), hf = 1u << (sh - 1u);
        if (rem > hf || (rem == hf && (hm & 1u) != 0u)) hm += 1u;
        return s | hm;
    }
    uint hm = m >> 13u, rem = m & 0x1FFFu;
    uint r = s | (uint(he) << 10u) | hm;
    if (rem > 0x1000u || (rem == 0x1000u && (hm & 1u) != 0u)) r += 1u;   // a carry rounds into the exponent
    return r;
}

#endif // MD_GPU_SLANG_INCLUDED
