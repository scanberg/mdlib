# md_gpu audit — Metal backend

Scope: `src/core/md_gpu.h`, `src/core/md_gpu_metal.m`, `src/core/md_gpu_builtin_msl.inl`,
`src/shaders/md_gpu.slang`, `unittest/test_gpu.c`, `unittest/shaders/gpu_test.slang`,
with `src/core/md_gpu_vulkan.c` read where needed to tell "Metal-only" from "both".

Nothing here was executed. The shell available to me is a Linux VM with no GPU, so
every claim below is either read off the source or off build artefacts in
`build/gen/`. Where the evidence is an artefact rather than the source, that is
called out, because `build/gen/` was generated 2026-08-10 and the shaders were
last touched 2026-08-17.

---

## 0. The suite has never been built in this tree

`CMakeLists.txt:479` wires `gpu_test.slang` through `compile_gpu_shaders`, and
`unittest/CMakeLists.txt:51` adds `test_gpu.c`. But `build/gen/` contains
`topo_*.metal` and `gto_*.metal` and no `gpu_test_*` at all — not the `.metal`,
not the `.metallib`, not `gpu_test_shaders.inl`. The GPU tests were added after
the last configure and have not been compiled, let alone run.

That matters for reading the rest of this document: none of the existing tests
have ever passed on Metal, so "the tests pass" is not evidence of anything yet.

---

## 1. Blocking correctness

### 1.1 The fixture never declared `group_size`, so every Metal kernel ran 1 thread per group — FIXED

`md_gpu_kernel_desc_t::group_size` is documented as *required on Metal*, and
`md_gpu_kernel_create` silently defaults it to `{1,1,1}` when left zero
(`md_gpu_metal.m:1498`). The `GPU_KERNEL` fixture macro never set it.

Consequence on Metal: `md_gpu_grid_1d(N, 64)` dispatches `N/64` threadgroups of
**one** thread, so `fill` writes 1/64 of the buffer and every data-carrying test
fails. `kernel_info_recovers_group_size` fails outright. On Vulkan the group size
is recovered from the SPIR-V, so the same source passes — which is exactly the
kind of divergence that makes this hard to spot.

Fixed in `unittest/test_gpu.c`: `GPU_KERNEL` now takes `(gx, gy, gz)` and each
kernel declares the size matching its `[numthreads]`.

**Recommendation:** make the zero case an error on Metal rather than a silent
`{1,1,1}`. A kernel that dispatches 1/64 of its threads produces plausible-looking
partial output, which is the worst failure mode available.

### 1.2 RETRACTED — bindless textures are fine on Metal

An earlier revision of this document claimed Slang hoists `DescriptorHandle<T>`
out of the argument struct into a separate `[[texture(0)]]` binding that md_gpu
never sets. **That was wrong.** The evidence was `build/ext/mdlib/gen/*.metal`
files dated 2026-08-10, generated from the *pre-migration* shaders that still used
`ConstantBuffer<RootArgs>` directly. The freshly generated `*_main.metal` files
have no texture binding at all, and `gpu_test_tex_write.metal` shows the handle
staying exactly where the design says it does:

```
struct TexArgs_0 {
    uint dim_x_0; uint dim_y_0; uint dim_z_0; float scale_0;
    texture3d<float, access::read_write> tex_0;      // in the struct, 8 bytes
};
```

So the design works: the handle is an `MTLResourceID` inside a Metal 3 argument
buffer, no member shifts, and no `setTexture:` is needed. `md_gpu_tex_sampled`
returning the same handle is correct too.

What actually broke the texture tests is §1.2a and §1.2b below — they broke
*every* kernel, textures included.

### 1.2a The root buffer held the argument struct, not a pointer to it — FIXED

This is the one that explains the failures. Slang does **not** collapse the
indirection on Metal. `gpu_test_fill.metal`:

```
struct FillRoot_0 { FillArgs_0 device* args_0; };
[[kernel]] void fill(uint3 tid [[thread_position_in_grid]],
                     FillRoot_0 constant* fill_root_0 [[buffer(0)]])
{
    FillArgs_0 device* _S1 = fill_root_0->args_0;    // <- dereferences buffer 0
    ...
}
```

The bound buffer must hold an **8-byte device address**. `md_mtl_encode_dispatch`
bound the argument struct itself, so the shader read `FillArgs.{n, base}` — for
the failing test, `{128, 5000}` — as a pointer and wrote through it. Hence
`5049 vs 0`.

`md_mtl_place_args` had computed the correct address into `op.arg_addr` all
along; nothing ever bound it. That dead field was the tell, and I read it the
wrong way round the first time.

Fixed: `md_mtl_place_args` now allocates an 8-byte **root cell** alongside the
argument block, writes the block's address into it, and
`md_mtl_encode_dispatch` binds that cell with an ordinary `setBuffer:`. The
address is stable across a graph replay because `md_gpu_graph_args` patches the
block in place rather than moving it, and for a captured launch both the block
and its cell come from graph-owned pages.

The first attempt used `setBytes:` instead, on the reasoning that an 8-byte
inline constant is exactly what that call is for. It is not reliable here:
called repeatedly at one index inside a single encoder, the second call did not
take effect, and the second of two back-to-back `fill` launches ran with the
first one's arguments — the upper half of `pointer_arithmetic_subranges` stayed
zero. An explicit cell has no shared mutable state: distinct offset, distinct
bytes, lifetime tied to the page like everything else.

Related, and removed for the same reason: a `bound_pso` cache that skipped
`setComputePipelineState:` when consecutive launches used the same kernel. It
silenced a debug-layer perf hint and bought a few CPU cycles, against a failure
mode of "this dispatch silently uses the previous launch's state". Every
dispatch now rebinds pipeline and root cell unconditionally.

Two consequences handled with it:

- **The built-in `make_grid` MSL took the struct by reference** — it was the one
  kernel the old binding was correct for, which is why indirect launches worked.
  Rewritten to take `constant MdGpuRoot&`.
- **Argument pages must now be resident.** The struct is reached by raw
  `gpuAddress` with no `setBuffer:` to imply residency. Arena and graph pages
  went through `md_mtl_make_resident`, but on the pre-macOS-15
  `useResources:` fallback `dev->live_bufs` tracked only pool blocks.
  `md_mtl_page_create` / `_destroy` now maintain it too.

### 1.2b The root buffer index is not always 0 — FIXED

`#define MD_MTL_ARG_BUFFER_INDEX 0` was hardcoded. Slang assigns Metal buffer
indices **per file, in declaration order of the push-constant globals** — not per
entry point. Every production shader is a single-entry file and therefore always
got 0, which is why this never surfaced. `gpu_test.slang` declares eight roots:

```
fill 0   scale_add 1   sum_reduce 2   tex_write 3
tex_read 4   layout_probe 5   tex_probe 6   bump 7
```

So `fill` had the right slot with the wrong contents, and everything else had
nothing bound at all. `gpu_volume.slang` is affected the same way
(`eval_field 0`, `compact_above 1`, `gather 2`).

Fixed: `md_gpu_kernel_create` now builds the pipeline with
`MTLPipelineOptionBindingInfo` and reads the index out of
`MTLComputePipelineReflection.bindings`, storing it per kernel. The constant
survives only as the fallback.

This one is worth a note in `md_gpu.h`: the index being a per-file property means
a shader author can silently change another kernel's binding by adding a root
global above it. Reflection makes that harmless, but it is surprising.

### 1.2c Compute encoders were serial, but the backend barriers them — FIXED

After §1.2a and §1.2b, 38 of 47 tests passed. The remaining 9 all had one of two
shapes, and one cause.

```
static void md_mtl_barrier_if_needed(md_gpu_stream_t s) {
    ...
    [s->compute_enc memoryBarrierWithScope:MTLBarrierScopeBuffers | MTLBarrierScopeTextures];
}
```

`[cmd computeCommandEncoder]` returns an `MTLDispatchTypeSerial` encoder.
`memoryBarrierWithScope:` is only defined on `MTLDispatchTypeConcurrent`; Apple
documents calling it on a serial encoder as a programmer error. It does not
fail — it corrupts execution ordering.

Every failure fits:

| Shape | Tests |
|---|---|
| A blit immediately before a compute dispatch — the barrier lands on a freshly opened encoder and the dispatch is dropped, writing nothing | `arg_struct_layout_matches_shader` (all `0xEEEEEEEE` sentinel), `texture_handle_reaches_shader`, `texture_host_write_then_kernel_read`, `make_grid_and_indirect_launch` (`sum_reduce` after `memset`), `host_callback_observes_preceding_readback`, `concurrent_streams_from_threads` |
| A long chain of dependent dispatches — barriers ignored, dispatches overlap, lost updates | `program_order_chain` (40 steps, 6 landed), `graph_replayed_many_times` (16 replays, 5 landed), `program_order_across_submissions` |

And every kernel that passed was the *first* op in its command buffer, where
`md_mtl_stream_ensure_cmd` clears `needs_barrier` so no barrier is emitted:
`fill`, `tex_write`. `tex_write` passing while `tex_read` failed had nothing to
do with textures — `tex_read` simply has a `memcpy_to_tex` blit in front of it.

**Fix, second attempt.** The first attempt switched the encoder to
`MTLDispatchTypeConcurrent` to make the barrier legal. That was the wrong half:
it regressed `pointer_arithmetic_subranges` — two independent `fill` launches —
to writing nothing at all, and it makes correctness depend on md_gpu placing a
barrier before every dependent dispatch, by hand, forever.

The encoder stays **serial** and the barrier call is **deleted**. That is what
md_gpu.h already promises:

> Launches into one stream are strictly ordered and never overlap, exactly as in
> CUDA — there is no flag to opt out.

A serial encoder orders its dispatches, Metal orders encoders within a command
buffer, and a queue runs command buffers in commit order. Every ordering the
header guarantees holds structurally, so the Metal backend now emits no barriers
at all. `md_mtl_barrier_if_needed` survives as an empty function only to keep the
call sites symmetric with `md_gpu_vulkan.c`, where the barriers are real.

Consequence, and it is a real one: `md_gpu_stream_set_auto_order(false)` and
`md_gpu_stream_barrier` cannot do anything on Metal — ordering is a property of
the encoder, not something that can be switched off per region. Both are no-ops.
That errs safe: code written against the escape hatch stays correct, it just does
not get faster. Making them real means concurrent encoders and hand-placed
barriers, which is the arrangement that silently mis-ordered this backend for its
entire existence. The header already says "Expect never to use this" — it should
also say the hatch is Vulkan-only.

### 1.2d `enable_validation` was accepted and ignored — FIXED

`md_gpu_device_desc_t::enable_validation` appeared nowhere in
`md_gpu_metal.m`. The unit-test fixture sets it. Metal API validation reports a
serial-encoder `memoryBarrierWithScope:` immediately, so §1.2c would have been a
one-line diagnostic instead of nine silent failures.

It cannot be enabled from inside the process — Metal reads
`METAL_DEVICE_WRAPPER_TYPE` before the first Metal call. `md_gpu_device_create`
now logs how to turn it on rather than letting the flag look effective. Worth
running the suite as:

```
METAL_DEVICE_WRAPPER_TYPE=1 MTL_DEBUG_LAYER=1 ./build/bin/md_unittest --filter='gpu.*'
```

### 1.2e Nothing ordered memory across command buffers — FIXED

With 1.2a–d in, `program_order_chain` (40 dependent dispatches, one encoder)
passed while `program_order_across_submissions` (the same chain split across 8
submissions) lost two of its eight steps — `1 << 6` where `1 << 8` was expected.
That isolates the failure to the command-buffer boundary, where
`md_mtl_stream_ensure_cmd` asserted:

> Metal executes command buffers on a queue in the order they were committed,
> and each stream owns its queue, so work in a new command buffer is already
> ordered after the previous one.

Metal does start them in commit order. But the *memory* ordering between command
buffers is driven by hazard tracking, and this backend gives Metal nothing to
track — buffers are reached by raw `gpuAddress` and declared only through a
device-wide `MTLResidencySet`, which handles residency and explicitly not
hazards. Metal cannot know command buffer N+1 reads what N wrote.

This is the price of the design's best property. "No per-dispatch resource
declarations" is exactly what makes the API pleasant, and it is also what
switches off every implicit barrier Metal would otherwise insert. Inside a
command buffer that costs nothing, because serial dispatch and encoder
boundaries order things structurally. Across submissions, the dependency has to
be stated.

Fixed: a new command buffer now waits on the stream's own previous signal value.

```objc
if (s->submitted_value > 0) {
    [cb encodeWaitForEvent:s->timeline value:s->submitted_value];
}
```

One event wait per command buffer, only at a submission boundary. This is the
header's central promise — "work issued into a stream executes in issue order,
and every operation observes all writes made by the operations before it in that
stream" — so it is not optional. Concurrency still comes from other streams,
which are other queues.

Both the file header and the residency note now say this outright, because it is
non-obvious and the next person to add a code path here needs to know that
implicit sync is unavailable on principle rather than by accident.

### 1.2f Patched graph arguments were invisible to replay — FIXED

`graph_capture_and_replay` captures `fill` + `scale_add(mul=1, add=7)`, replays,
patches the second launch's block through `md_gpu_graph_args` to `mul=2,
add=1000`, and replays again. The second replay produced `i + 7` — the original
arguments.

The patch itself was fine: `EXPECT_EQ(7u, patched->add)` passed, so
`md_gpu_graph_args` returned the right memory, and the CPU wrote to it. And the
*first* replay was correct, so CPU→GPU visibility works in general. What differs
the second time is that the GPU had already read that page in an earlier command
buffer.

This is 1.2e again in another costume. Metal invalidates GPU-side caches for
resources it tracks; this backend declares nothing per dispatch, so a CPU write
to a page the GPU has previously read has nothing to trigger the invalidate. The
result is the worst possible shape: correct once, then silently stale forever.

Fixed by not binding the shadow at all. Capture now records only the argument
*bytes*, and `md_gpu_graph_launch` copies them into fresh arena memory —
allocated exactly as a non-captured launch would — before each dispatch. Memory
the GPU has never seen cannot be stale, so the whole question is moot rather than
carefully reasoned about.

The cost is a memcpy of a small struct per dispatch per replay. Nothing a graph
actually saves is given up: no re-recording, no pipeline lookup, no encoder
rebuild. The file header's claim that argument blocks "are never rebuilt" is now
wrong and should be reworded — they are rebuilt, cheaply, and that is the point.

A knock-on: the graph's argument pages are now pure CPU shadow. They are still
allocated as `MTLBuffer`s, which is harmless but no longer necessary; a plain
host allocation would do.

### 1.2g Encoders inside a command buffer were never ordered — FIXED

Reported as artifacts in the GTO density and MO volumes: mostly correct, with
scattered blocks wrong. Block-shaped is the tell — the GTO kernels are
`[numthreads(8,8,8)]` and screen per threadgroup, so a whole 8³ block is the
unit that goes bad.

Not the kernels. `git diff` on `eval_gto_density.slang` / `eval_gto_mo.slang`
across the refactor is pure argument plumbing — the `MdGpuRoot` migration and
moving the texture into the struct. No maths changed. And `md_gto.c` declares
`group_size` 8,8,8 matching `[numthreads]`, with the grid divided by 8 in each
axis, so the launch geometry is right too.

The hazard is in `veloxchem.cpp`'s readback path, and it is generic:

```c
md_gto_gpu_density_launch(gpu_stream, &desc);   // compute encoder, writes gpu_volume
...
md_gpu_memcpy_from_tex_async(rb, gpu_volume, &region, size, gpu_stream);  // blit encoder, reads it
md_gpu_stream_sync(gpu_stream);
```

Encoders inside one command buffer are **not** guaranteed to run one after
another. Metal may overlap them, and decides whether it may from hazard
tracking. The texture arrives through a bindless handle inside the argument
struct, so nothing declares it to either encoder — as far as Metal can tell the
dispatch and the blit are independent, and the blit may start copying voxels the
dispatch has not written yet.

Fixed with an `MTLFence` per stream: the encoder being closed calls
`updateFence:`, the next one calls `waitForFence:`. That is exactly what
`MTLFence` is for — ordering untracked resource access across encoders — and it
is far cheaper than splitting into separate command buffers. The fence is scoped
to a command buffer, so `fence_valid` is re-armed in `md_mtl_stream_ensure_cmd`
and no wait is issued before the first update.

Covered by the hazard matrix in §4.

### 1.3 A captured host→device copy records a pointer into recycled scratch

Argument blocks get their own graph-owned pages (`md_mtl_place_args`, the
`s->capture` branch). Staging does not:

- `md_gpu_memcpy_async` host→device during capture — `md_gpu_metal.m:1044`
- `md_gpu_upload_begin` / `md_gpu_upload_end` during capture — `md_gpu_metal.m:1115`

Both call `md_mtl_arena_alloc`, which draws from the **stream's** transient arena.
`md_mtl_arena_retire` marks those pages reusable once the capturing stream passes
the submission, and a graph is by definition replayed later. So the recorded blit
reads whatever the arena has since been used for.

This is silent — the copy still happens, with the wrong bytes.

**Fix:** during capture, stage into a graph-owned page (`g->arg_pages` already has
the right lifetime), or reject host-source copies inside a capture.

### 1.4 The device→host fast path does not check `s->capture`

```c
/* md_gpu_metal.m:1052 */
if (sblk->host && !s->has_work && md_mtl_stream_completed(s) >= s->submitted_value) {
    memcpy(dst, (uint8_t*)sblk->host + src_off, size);
    return true;
}
```

The host→device path four lines above correctly guards with `!s->capture`. This
one does not, and `md_mtl_did_op` deliberately leaves `has_work` false during
capture — so with a `MD_GPU_MEM_HOST_READ` source the copy executes *at capture
time* and is not recorded at all. Replaying the graph never performs it.

### 1.5 `md_gpu_memcpy_from_tex_async` into host memory is dropped when captured

`md_gpu_metal.m:1413` records the texture→buffer blit, then returns. The delivery
callback that copies staging into the caller's host buffer is only registered on
the non-capture path (`:1433`). A graph containing a texture readback therefore
never writes the host buffer, with no error.

### 1.6 Readback staging can be recycled before it is delivered

The device→host path stages into the arena, records a sync, and registers
`md_mtl_d2h_finish` against that sync value. `md_mtl_arena_retire` stamps the same
page with the *same* value. So the moment the GPU passes it, the page is eligible
for reuse by `md_mtl_arena_alloc` **and** the callback is eligible to run — and
nothing orders those two. A subsequent upload on the same stream, before the next
`md_gpu_device_poll`, overwrites the bytes the callback is about to copy out.

`md_gpu_stream_sync` masks it by calling `md_mtl_poll_internal(dev, false)`, which
is why the existing suite doesn't see it. New test:
`readback_survives_subsequent_staged_work` interleaves an upload between the
readback and the poll.

**Fix:** retire staging pages at `sync.value + 1`, or hold a reference on the page
in `md_mtl_d2h_t` and release it in the callback.

### 1.7 The allocation registry is read without its lock — both backends

`device_mutex` is documented as covering "pools, registry, textures, retire lists"
and `md_gpu_malloc` / `md_gpu_free` take it. These callers do not:

- Metal: `md_gpu_memcpy_async:1029-1030`, `md_gpu_memset_async:1079`,
  `md_gpu_upload_begin:1106`, `md_gpu_upload_end:1131`,
  `md_gpu_launch_indirect:1630`, `md_gpu_memcpy_to_tex_async:1360`,
  `md_gpu_memcpy_from_tex_async:1403`
- Vulkan: `md_gpu_vulkan.c:1844, 1845, 1914, 1957, 1985`

`md_mtl_registry_insert` `memmove`s the array under the lock. A concurrent
`md_gpu_memcpy_async` binary-searching it is a straight data race, and the header
explicitly promises "different streams may be used concurrently from different
threads" and "md_gpu_malloc / md_gpu_free ... are thread-safe". New test:
`concurrent_streams_from_threads` (each thread owns its stream and its memory, so
any failure is md_gpu's shared state). Run it under TSan.

Related, Metal-only: `md_mtl_devices[]` / `md_mtl_device_count` are plain globals
mutated by `md_gpu_device_create` / `_destroy` with no lock, and read lock-free by
`md_mtl_find_anywhere` and `md_mtl_tex_device`.

### 1.8 `md_gpu_pool_destroy` does not wait for work in flight

It destroys every block unconditionally, `in_use` or not
(`md_gpu_metal.m:873-874`). Buffers reached through a 64-bit `gpuAddress` are not
retained by the command buffer that uses them — that is the whole point of the
residency set — so releasing one with a dispatch in flight is a use-after-free.

`md_gpu_stream_destroy` and `md_gpu_device_destroy` both sync first. `pool_destroy`
is the odd one out, and the header says nothing either way. New test:
`pool_destroy_with_work_in_flight`.

**Fix:** either sync the streams that touched the pool, or document that the
caller must. Given the rest of the API's stream-ordered-and-safe posture, syncing
is the consistent choice.

### 1.9 `md_gpu_device_destroy` does not own kernels or graphs

The header: *"Waits for all streams to go idle, then destroys the device and
everything created from it."*

`md_gpu_device` tracks pools, streams, textures and samplers. It does not track
kernels or graphs. So an undestroyed kernel leaks its `MTLComputePipelineState`,
and — worse — `md_gpu_kernel_destroy` called after `md_gpu_device_destroy`
dereferences `k->device->alloc`, which is freed memory. Same for
`md_gpu_graph_destroy`. New test:
`device_destroy_reclaims_undestroyed_objects`.

### 1.10 Fixed-size arrays truncate silently

| Site | Limit | On overflow |
|---|---|---|
| `md_mtl_retire_t.waits[8]` | 8 streams | texture released while a 9th stream may still use it — **silent** |
| `MD_MTL_MAX_PENDING_WAITS` | 16 | `md_mtl_fail`, then the wait is *dropped* and `md_gpu_stream_wait` returns void — the caller cannot tell |
| `md_mtl_timeline_snapshot_t[16]` | 16 streams | falls back to a fresh read, reintroducing the exact ordering bug the snapshot exists to prevent |

New test: `wide_cross_stream_fan_in` (20 producer streams into one consumer).

### 1.11 Queued cross-stream waits can be lost

`md_gpu_stream_wait` submits only `if (s->has_work)`, then queues the wait for
`md_mtl_stream_ensure_cmd` to encode at the head of the *next* command buffer. But
`md_mtl_stream_submit` returns early when `!s->has_work` **without releasing
`s->cmd`**, and `ensure_cmd` returns an existing `s->cmd` untouched. So if a
command buffer exists with no work — reachable through `md_gpu_capture_begin`,
which calls `md_mtl_stream_submit` unconditionally — the queued waits are never
encoded and the cross-stream dependency silently disappears.

### 1.12 A stream's timeline can dangle

`s->wait_events[i]` is an `__unsafe_unretained id<MTLSharedEvent>` belonging to
another stream. `md_gpu_stream_destroy` releases that event
(`md_mtl_stream_destroy_internal:745`) with no check for waiters. Destroy a
producer stream while a consumer holds a queued wait on it and the consumer
encodes against a released object.

---

## 2. Design

### 2.1 The pool is a block cache, not a suballocator

`md_gpu_pool_desc_t::block_size` is documented as a "suballocation granularity
hint" and is stored and never read. Every `md_gpu_malloc` is one whole `MTLBuffer`
rounded up by `md_mtl_next_pow2`. Two consequences:

- **Waste.** A 1.4 GB volume becomes a 2 GB buffer. Volume grids are rarely powers
  of two, so this is the common case, not the corner.
- **Cost per allocation.** Each one is a `newBufferWithLength:` plus a residency
  `addAllocation:`+`commit`, plus an O(n) insert into the sorted registry. Small
  allocations are expensive in a way `bytes_in_use` does not show — it reports the
  rounded capacity, not the requested size, so a 5 KB request reads as 8 KB.

Either implement suballocation or drop `block_size` and say plainly that a pool
caches whole device allocations.

### 2.2 Residency is committed once per allocation

`md_mtl_make_resident` / `md_mtl_end_residency` call `commit` on every single
allocation and free. `MTLResidencySet` commit is the expensive operation in that
API; it is meant to be batched. Accumulate additions and commit once, lazily,
before the next submit.

### 2.3 Both wait paths spin

```c
/* md_gpu_metal.m:720 */
/* MTLSharedEvent has no blocking wait in C; spin with a short yield. */
while (...) md_thread_sleep(0);
```

The comment is not correct. `MTLSharedEvent` has
`-waitUntilSignaledValue:timeoutMS:` (macOS 12+), and `MTLCommandBuffer` has
`-waitUntilCompleted` and `-addCompletedHandler:`. `md_gpu_stream_sync`,
`md_gpu_sync_wait` and the loop in `md_gpu_device_destroy` all burn a core for the
duration of the GPU work. On a battery-powered laptop running a long volume
evaluation this is very visible.

### 2.4 Graphs pin nothing

A graph records raw `id<MTLBuffer>` / `id<MTLTexture>` (`md_mtl_op_t`,
`__unsafe_unretained`) with no ownership and no liveness check. Free the
allocation, or `md_gpu_pool_reset` the pool it came from, and replaying the graph
is a use-after-free. Nothing in the header's Graphs section warns about this,
while the same section actively encourages long-lived graphs ("the graph is never
re-recorded").

At minimum document it. Better: have the graph retain the blocks it references, or
record `(pool, offset)` and re-resolve at launch.

### 2.5 `md_gpu_stream_record` does not do what the header says

> *"Returns the none sync if nothing has been issued since the last record."*

Both backends return `s->submitted_value`, i.e. the last submission, so a second
`record` with nothing in between returns the same valid sync rather than none.
The implementations agree with each other, so this is a documentation bug — but
`md_gpu_launch_host_fn` builds on `record`, and "fire after the previous batch"
versus "fire immediately" is a meaningful difference for a caller reading the
header.

### 2.6 `md_gpu_capture_next_index` counts the wrong thing

> *"The ordinal the next md_gpu_launch into this stream will be given."*

It returns `capture->launches.count`, and `md_mtl_place_args` only pushes a launch
record when `args_size > 0`. A zero-argument launch shifts every subsequent index,
so `md_gpu_graph_args(g, i)` patches the wrong block.

### 2.7 No validation on texture copies

`md_gpu_memcpy_to_tex_async` / `_from_tex_async` never check `size` against the
region's volume — a short buffer is an out-of-bounds read with no diagnostic. And
`md_mtl_resolve_region` computes `w - r->offset[0]` unsigned, so an offset past
the extent underflows to a ~4-billion-texel extent. Both are cheap to check. The
new `texture_region_subrange_roundtrip` test exercises the happy path; the
validation is still missing.

### 2.8 Smaller things

- **The built-in `make_grid` MSL broke the header's own rule — FIXED.** It
  declared `uint3 local` against a C mirror using `uint32_t local[3]`
  (`md_mtl_make_grid_args_t`), and the two structs were 32 and 48 bytes. It
  worked only because every member the kernel reads happens to fall inside the
  first 32 bytes. Now `md_gpu_uint3` on the C side and no trailing pad on the
  MSL side; both are 32 bytes with matching offsets, verified by compiling the
  mirror under each backend define.
- **`check_gpu_arg_layout.py` did not catch either §1.2a or §1.2b.** It compares
  *argument struct* layouts across targets, which were fine. Neither the root
  indirection nor the buffer index is a layout property, so nothing in the check
  models them. Worth extending: it already parses the MSL, so asserting that the
  root struct's sole member is a pointer, and recording the assigned
  `[[buffer(N)]]`, is cheap and would have caught both at build time.
- **Dead device info.** `timestamp_period_ns_num/den` are hardcoded 1/1 and there
  is no timestamp API to use them with. `preferred_group_multiple` is hardcoded 32
  rather than read from a pipeline's `threadExecutionWidth`.
- **`bytes_cached` overstates.** It is `reserved - in_use`, but a block freed on a
  still-running stream is not reusable by a different stream. The header says
  "reusable without a new device allocation".
- **Leaks on failure paths.** `md_gpu_kernel_create` leaks the retained `fn`/`pso`
  if the following `md_alloc` fails; `md_mtl_stream_create_internal` leaks the
  command queue if `newSharedEvent` fails.

---

## 3. Test coverage

The existing 21 tests are well chosen — program order across submissions, the
callback-ordering race, stream-ordered pool reset. What they did not reach:

| Gap | New test |
|---|---|
| `md_gpu_uint3` / `uint2` / `float4` mirrors — the header's largest section, untested | `arg_struct_layout_matches_shader` |
| Whether a `DescriptorHandle` reaches the shader at all, separately from the round trip | `texture_handle_reaches_shader` |
| Staging larger than one arena page | `upload_larger_than_one_arena_page` |
| Readback staging vs. subsequent uploads (1.6) | `readback_survives_subsequent_staged_work` |
| Copies and fills inside a graph, replayed after the arena has churned (1.3) | `graph_replays_captured_copies_and_fills` |
| Texture regions — origin, extent, row/slice pitch | `texture_region_subrange_roundtrip` |
| Samplers | `sampler_create_and_destroy` |
| `md_gpu_tex_sampled` on a dual-usage texture | `texture_storage_and_sampled_handles` |
| More streams than the fixed wait/snapshot arrays (1.10) | `wide_cross_stream_fan_in` |
| The documented threading contract (1.7) | `concurrent_streams_from_threads` |
| Pool destroy with work in flight (1.8) | `pool_destroy_with_work_in_flight` |
| Device destroy as the owner of last resort (1.9) | `device_destroy_reclaims_undestroyed_objects` |
| NULL / zero-size / empty-grid tolerance | `null_and_zero_arguments_are_tolerated` |
| `md_gpu_host_ptr` and pointer arithmetic on host-visible memory | `host_pointer_tracks_device_pointer` |
| `md_gpu_last_error` after a failure; double `upload_begin` | `last_error_describes_the_failure` |
| Device info sanity; kernel group size within the device limit | `device_info_is_sane` |
| Indirect grid derived from a zero count | `indirect_launch_with_zero_count` |
| Pool reuse handed to an unordered stream | `pool_reuse_across_streams_waits_for_completion` |

Two new shader entries back these: `layout_probe` and `tex_probe` in
`unittest/shaders/gpu_test.slang`, added to the `ENTRIES` list at
`CMakeLists.txt:479`.

### Verification done without a GPU

- `unittest/test_gpu.c` passes `gcc -fsyntax-only -std=c11` against the real
  headers with a stubbed `gpu_test_shaders.inl`.
- The two new C mirrors were checked field by field against the SPIR-V-scalar and
  MSL layout rules, compiled once with `MD_GPU_BACKEND_METAL` and once with
  `MD_GPU_BACKEND_VULKAN`. All offsets and both `sizeof`s agree:

  ```
  METAL   dim@0  scale@16  pair@24  v4@32  dst@48  sizeof 64
  VULKAN  dim@0  scale@12  pair@16  v4@24  dst@40  sizeof 48
  tex_probe (both):  n@0  dst@16  tex@24  marker@32  sizeof 40
  ```

---

## 4. The hazard matrix

Four of the bugs above (§1.2c, §1.2e, §1.2f, §1.2g) are one bug wearing different
hats: **md_gpu.h promises that every operation observes the writes of the
operations before it in the same stream, and neither backend gets that for
free.** Metal decides whether adjacent encoders may overlap from hazard tracking
it is never given; Vulkan does whatever `vkCmdPipelineBarrier` the backend emits,
and a missing or over-narrow stage/access mask fails in exactly the same places.

The whole class went undetected through a passing suite and surfaced as visible
artifacts in a GTO volume. That is a test-design failure, not bad luck, and it
came down to two things:

- **Size.** A producer has to run long enough for an unordered consumer to get
  ahead of it. The original texture tests use an 8³ volume — one threadgroup,
  finished microseconds before a blit could possibly overtake it. The hazard was
  always present; nothing ran long enough to expose it. The matrix uses 128³
  volumes and 1 M-element buffers.
- **Rounds.** Each round writes a different value, so a stale read reports the
  *previous round's* value. "Expected 3000, got 2000" says ordering; "expected
  3000, got 0" says something else entirely. A single-round test cannot tell
  those apart.

Nothing in these tests mentions fences, barriers or encoders — they are written
against md_gpu.h's guarantee, so they constrain both backends equally.

| Test | Producer → consumer | Separation exercised |
|---|---|---|
| `hazard_dispatch_to_blit_buffer` | dispatch → transfer, buffer | encoder boundary |
| `hazard_dispatch_to_blit_texture` | dispatch → transfer, texture | encoder boundary — the GTO readback, reduced |
| `hazard_blit_to_dispatch_buffer` | transfer → dispatch, buffer | encoder boundary, other direction |
| `hazard_blit_to_dispatch_texture` | transfer → dispatch, texture | encoder boundary, other direction |
| `hazard_memset_to_dispatch` | fill → dispatch | fill is a transfer op too |
| `hazard_dispatch_to_dispatch_via_texture` | dispatch → dispatch, texture | within one encoder (serial dispatch on Metal, a barrier on Vulkan) |
| `hazard_alternating_encoder_chain` | 12× (transfer, dispatch) | many transitions in one command buffer; one missed link shows at the end |
| `hazard_across_submission_boundary` | dispatch → transfer, split by `flush` | submission boundary (§1.2e) |
| `hazard_across_streams_via_sync` | dispatch → transfer, two streams | `record` / `wait` carries a texture dependency |
| `hazard_inside_a_replayed_graph` | dispatch → transfer, captured | replay re-encodes, so the ordering must be re-inserted |

The last one is worth keeping: replay walks a different code path from recording,
and inserting the ordering once during recording and forgetting it during replay
is an easy and entirely silent mistake.

Two gaps the matrix does not close, both cheap to add later if they ever bite:
a dispatch writing a texture that a *sampled* read consumes (only storage access
is covered), and hazards on the transfer-kind stream where a backend may pick a
dedicated DMA queue.

### Suggested order on the machine

```
cmake -S ext/mdlib -B ext/mdlib/build -DMD_ENABLE_GPU=ON -DCMAKE_BUILD_TYPE=Debug
cmake --build ext/mdlib/build --target md_unittest -j
./ext/mdlib/build/bin/md_unittest --filter='gpu.*'
```

Read the failures in this order, because the earlier ones invalidate the later:

1. `launch_writes_buffer` — the simplest kernel. Confirms §1.1/§1.2a/§1.2b.
2. `arg_struct_layout_matches_shader` — if this fails nothing else means anything.
3. `texture_handle_reaches_shader` — `dst[0]` reports whether the struct survived
   the handle, the texture readback whether the handle resolved; they fail
   independently.
4. `make_grid_and_indirect_launch` — the rewritten built-in.
5. everything else.

Then re-run under `-fsanitize=thread` for `concurrent_streams_from_threads` and
`-fsanitize=address` for the lifetime tests.
