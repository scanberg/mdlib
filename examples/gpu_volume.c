/*
gpu_volume.c -- a worked end-to-end example of md_gpu.

It runs the shape viamd actually runs, in miniature:

    evaluate a scalar field into a 3D volume
      -> compact the voxels above a threshold, counting on the device
      -> process exactly that many with an indirect launch
      -> read the results back without blocking the calling thread

and then re-issues the compaction with different thresholds, which is what a
slider drag does.

Along the way it touches most of the API: streams, stream-ordered allocation,
zero-copy uploads, GPU addresses and address arithmetic, textures and their
shader handles, cross-stream synchronisation, indirect launches driven by a
device-side count, explicit ordering with stage barriers, and host callbacks.

Build: part of mdlib when MD_ENABLE_GPU is on. Run it directly; it prints what
it did and returns non-zero if a check fails.
*/

#include <core/md_gpu.h>
#include <core/md_allocator.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gpu_volume_shaders.inl"

/* ---------------------------------------------------------------------------
   Argument structs, mirroring examples/shaders/gpu_volume.slang.
   The shader's volume dimension travels as a uint4 -- argument structs may not
   contain 3-vectors -- and md_gpu_uint4 carries the 16-byte alignment the ABI
   requires, so this one C struct is correct on every backend. Pointer fields
   are md_gpu_addr_t: an allocation's .gpu goes in with no cast.
   --------------------------------------------------------------------------- */

typedef struct {
    md_gpu_uint4         dim;
    float                freq;
    uint32_t             _pad;
    md_gpu_storage_tex_t vol;
} eval_args_t;

typedef struct {
    md_gpu_uint4         dim;
    float                threshold;
    uint32_t             _pad0;
    md_gpu_storage_tex_t vol;
    md_gpu_addr_t        count;
    md_gpu_addr_t        indices;
    uint32_t             capacity;
    uint32_t             _pad1;
} compact_args_t;

typedef struct {
    md_gpu_uint4         dim;
    md_gpu_storage_tex_t vol;
    md_gpu_addr_t        count;
    md_gpu_addr_t        indices;
    md_gpu_addr_t        values;
} gather_args_t;

#define CHECK(cond, msg)                                                   \
    do {                                                                    \
        if (!(cond)) {                                                      \
            fprintf(stderr, "FAILED: %s (%s)\n", msg,                       \
                    md_gpu_last_error() ? md_gpu_last_error() : "-");     \
            return 1;                                                       \
        }                                                                   \
    } while (0)

enum { DIM = 32, VOXELS = DIM * DIM * DIM };

/* Delivered by a host callback once the readback has actually landed. */
typedef struct {
    const uint32_t* count;       /* host-readable memory, written by the GPU */
    const float*    values;
    int             fired;
} readback_t;

static void on_readback_complete(void* user) {
    readback_t* rb = (readback_t*)user;
    rb->fired = 1;
    /* Everything the stream had issued before md_gpu_launch_host_fn has
       completed, so the host-readable copies are populated by now. */
    printf("  [callback] readback complete: %u voxels above threshold, "
           "first value %.4f\n", *rb->count, *rb->count ? rb->values[0] : 0.0f);
}

static uint32_t cpu_count_above(float threshold) {
    uint32_t n = 0;
    for (int z = 0; z < DIM; ++z)
    for (int y = 0; y < DIM; ++y)
    for (int x = 0; x < DIM; ++x) {
        float px = ((float)x + 0.5f) / DIM * 2.0f - 1.0f;
        float py = ((float)y + 0.5f) / DIM * 2.0f - 1.0f;
        float pz = ((float)z + 0.5f) / DIM * 2.0f - 1.0f;
        if (expf(-4.0f * (px*px + py*py + pz*pz)) > threshold) n++;
    }
    return n;
}

int main(void) {
    /* -----------------------------------------------------------------
       Device
       ----------------------------------------------------------------- */
    md_gpu_device_desc_t dd = {0};
    dd.enable_validation = true;
    dd.label             = "gpu_volume example";

    md_gpu_device_t dev = md_gpu_device_create(&dd);
    if (!dev) {
        printf("no GPU device available: %s\n",
               md_gpu_last_error() ? md_gpu_last_error() : "-");
        return 0;    /* not a failure: there may simply be no GPU here */
    }

    md_gpu_device_info_t info;
    md_gpu_device_info(dev, &info);
    printf("device: %s (%s, max %u threads/group)\n",
           info.name, info.is_discrete ? "discrete" : "unified memory",
           info.max_threads_per_group);

    /* Two streams. Work in one is ordered; the two are independent until
       joined by a sync point. */
    md_gpu_stream_t compute  = md_gpu_stream_default(dev, MD_GPU_STREAM_COMPUTE);
    md_gpu_stream_t transfer = md_gpu_stream_default(dev, MD_GPU_STREAM_TRANSFER);

    /* -----------------------------------------------------------------
       Kernels. The generated descriptors carry each kernel's group size and
       argument-struct size, so nothing here repeats [numthreads].
       ----------------------------------------------------------------- */
    md_gpu_kernel_desc_t kd;
    kd = md_shader_gpu_volume_eval_field_kernel();    md_gpu_kernel_t k_eval    = md_gpu_kernel_create(dev, &kd);
    kd = md_shader_gpu_volume_compact_above_kernel(); md_gpu_kernel_t k_compact = md_gpu_kernel_create(dev, &kd);
    kd = md_shader_gpu_volume_gather_kernel();        md_gpu_kernel_t k_gather  = md_gpu_kernel_create(dev, &kd);
    CHECK(k_eval && k_compact && k_gather, "kernel creation");

    md_gpu_kernel_info_t ki;
    md_gpu_kernel_info(k_eval, &ki);
    printf("eval_field threadgroup: %ux%ux%u, %u-byte arguments\n",
           ki.group_size[0], ki.group_size[1], ki.group_size[2], ki.args_size);

    /* -----------------------------------------------------------------
       Resources. A texture is an object; what goes into an argument struct
       is its storage handle. Creation is stream-ordered and never waits.
       ----------------------------------------------------------------- */
    md_gpu_texture_desc_t td = {0};
    td.type   = MD_GPU_TEX_3D;
    td.format = MD_GPU_FORMAT_R32_FLOAT;
    td.usage  = MD_GPU_TEX_STORAGE;
    td.width  = DIM; td.height = DIM; td.depth_or_layers = DIM;
    td.label  = "field";
    md_gpu_texture_t vol = md_gpu_texture_create(compute, &td);
    CHECK(vol != NULL, "texture creation");
    const md_gpu_storage_tex_t vol_h = md_gpu_texture_storage(vol, 0);

    /* Stream-ordered allocation, by memory kind. Freeing is legal at any
       point, even with work in flight; the memory is reused without a fence. */
    md_gpu_addr_t count   = md_gpu_malloc(compute, MD_GPU_MEM_DEVICE, sizeof(uint32_t)).gpu;
    md_gpu_addr_t grid    = md_gpu_malloc(compute, MD_GPU_MEM_DEVICE, 3 * sizeof(uint32_t)).gpu;
    md_gpu_addr_t indices = md_gpu_malloc(compute, MD_GPU_MEM_DEVICE, VOXELS * sizeof(uint32_t)).gpu;
    md_gpu_addr_t values  = md_gpu_malloc(compute, MD_GPU_MEM_DEVICE, VOXELS * sizeof(float)).gpu;
    /* Results land in host-readable memory, read through .cpu once complete. */
    md_gpu_mem_t  rb_count  = md_gpu_malloc(compute, MD_GPU_MEM_HOST_READ, sizeof(uint32_t));
    md_gpu_mem_t  rb_values = md_gpu_malloc(compute, MD_GPU_MEM_HOST_READ, VOXELS * sizeof(float));
    CHECK(count && grid && indices && values && rb_count.cpu && rb_values.cpu, "allocation");

    /* A zero-copy upload, for data you would otherwise pack into a scratch
       buffer and memcpy. Here it seeds the second half of `values` -- plain
       byte arithmetic on the address -- purely to show the shape. */
    const md_gpu_addr_t upper = values + (VOXELS / 2) * sizeof(float);
    float* seed = md_gpu_upload_begin(transfer, upper, 16 * sizeof(float));
    CHECK(seed != NULL, "upload_begin");
    for (int i = 0; i < 16; ++i) seed[i] = (float)i;
    CHECK(md_gpu_upload_end(transfer), "upload_end");

    /* Join the two streams: compute waits for the transfer to land. */
    md_gpu_stream_wait(compute, md_gpu_stream_record(transfer));

    /* -----------------------------------------------------------------
       The pipeline. No barriers, no usage flags, no fences: everything
       issued into `compute` runs in order and sees the previous writes.
       ----------------------------------------------------------------- */
    const float threshold = 0.5f;
    const md_gpu_uint4 dim = {DIM, DIM, DIM, 0};

    eval_args_t ea = {0};
    ea.dim = dim; ea.freq = 4.0f; ea.vol = vol_h;
    CHECK(MD_GPU_LAUNCH(compute, k_eval, md_gpu_grid_for(k_eval, DIM, DIM, DIM), ea), "eval_field");

    CHECK(md_gpu_memset(compute, count, 0, sizeof(uint32_t)), "memset count");

    compact_args_t ca = {0};
    ca.dim = dim; ca.threshold = threshold; ca.vol = vol_h;
    ca.count = count; ca.indices = indices; ca.capacity = VOXELS;
    CHECK(MD_GPU_LAUNCH(compute, k_compact, md_gpu_grid_for(k_compact, DIM, DIM, DIM), ca), "compact_above");

    /* Turn the device-side count into an indirect grid sized for `gather`.
       The count is never read back to decide how much work to launch. */
    CHECK(md_gpu_make_grid(compute, grid, count, k_gather), "make_grid");

    gather_args_t ga = {0};
    ga.dim = dim; ga.vol = vol_h; ga.count = count; ga.indices = indices; ga.values = values;
    CHECK(md_gpu_launch_indirect(compute, k_gather, grid, &ga, sizeof(ga)), "gather");

    /* -----------------------------------------------------------------
       Non-blocking readback. Nothing here waits.
       ----------------------------------------------------------------- */
    static readback_t rb;
    rb.count  = (const uint32_t*)rb_count.cpu;
    rb.values = (const float*)rb_values.cpu;

    CHECK(md_gpu_copy(compute, rb_count.gpu, count, sizeof(uint32_t)), "readback count");
    CHECK(md_gpu_copy(compute, rb_values.gpu, values, VOXELS * sizeof(float)), "readback values");
    CHECK(md_gpu_launch_host_fn(compute, on_readback_complete, &rb), "host fn");

    /* In a real frame loop this is the only call you make; the callback fires
       on this thread, so touching OpenGL or ImGui from it is legal. */
    while (!rb.fired) {
        md_gpu_device_poll(dev);
    }

    const uint32_t expect = cpu_count_above(threshold);
    const uint32_t got    = *rb.count;
    printf("voxels above %.2f: gpu=%u cpu=%u %s\n",
           threshold, got, expect, got == expect ? "(match)" : "(MISMATCH)");
    CHECK(got == expect, "count matches the analytic result");
    for (uint32_t i = 0; i < got; ++i) {
        CHECK(rb.values[i] > threshold, "every gathered value is above the threshold");
    }

    /* -----------------------------------------------------------------
       A slider drag: re-issue the compaction with new thresholds. Issuing
       costs one argument-struct copy per launch; nothing is re-recorded
       or rebuilt. This loop also shows EXPLICIT ordering -- the caller
       places coarse stage barriers, with no resource lists.
       ----------------------------------------------------------------- */
    const float thresholds[] = {0.25f, 0.75f};
    md_gpu_stream_set_ordering(compute, MD_GPU_ORDER_EXPLICIT);
    for (int t = 0; t < 2; ++t) {
        /* The previous readback copy read `count`; the memset overwrites it. */
        md_gpu_barrier(compute, MD_GPU_STAGE_TRANSFER, MD_GPU_STAGE_TRANSFER);
        CHECK(md_gpu_memset(compute, count, 0, sizeof(uint32_t)), "memset count");
        md_gpu_barrier(compute, MD_GPU_STAGE_TRANSFER, MD_GPU_STAGE_COMPUTE);

        ca.threshold = thresholds[t];
        CHECK(MD_GPU_LAUNCH(compute, k_compact, md_gpu_grid_for(k_compact, DIM, DIM, DIM), ca), "compact_above");
        md_gpu_barrier(compute, MD_GPU_STAGE_COMPUTE, MD_GPU_STAGE_TRANSFER);

        CHECK(md_gpu_copy(compute, rb_count.gpu, count, sizeof(uint32_t)), "readback count");
        md_gpu_stream_sync(compute);

        const uint32_t n    = *(const uint32_t*)rb_count.cpu;
        const uint32_t want = cpu_count_above(thresholds[t]);
        printf("threshold %.2f: gpu=%u cpu=%u %s\n",
               thresholds[t], n, want, n == want ? "(match)" : "(MISMATCH)");
        CHECK(n == want, "re-issued chain honours the new threshold");
    }
    md_gpu_stream_set_ordering(compute, MD_GPU_ORDER_IMPLICIT);

    /* -----------------------------------------------------------------
       Teardown. Freeing is stream-ordered and never blocks; destroying the
       texture is legal even with work still in flight.
       ----------------------------------------------------------------- */
    md_gpu_memory_stats_t st;
    md_gpu_memory_stats(dev, MD_GPU_MEM_DEVICE, &st);
    printf("device memory: %llu in use in %u allocations, %llu reserved in %u chunks, %llu in textures\n",
           (unsigned long long)st.bytes_in_use, st.allocations,
           (unsigned long long)st.bytes_reserved, st.chunks, (unsigned long long)st.bytes_textures);

    md_gpu_free(compute, count);
    md_gpu_free(compute, grid);
    md_gpu_free(compute, indices);
    md_gpu_free(compute, values);
    md_gpu_free(compute, rb_count.gpu);
    md_gpu_free(compute, rb_values.gpu);
    md_gpu_texture_destroy(vol);
    md_gpu_kernel_destroy(k_eval);
    md_gpu_kernel_destroy(k_compact);
    md_gpu_kernel_destroy(k_gather);
    md_gpu_device_destroy(dev);      /* waits for idle internally */

    printf("OK\n");
    return 0;
}
