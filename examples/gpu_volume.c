/*
gpu_volume.c -- a worked end-to-end example of md_gpu.

It runs the shape viamd actually runs, in miniature:

    evaluate a scalar field into a 3D volume
      -> compact the voxels above a threshold, counting on the device
      -> process exactly that many with an indirect launch
      -> read the results back without blocking the calling thread

and then captures the same chain into a graph and replays it with a different
threshold, which is what a slider drag would do.

Along the way it touches every part of the API: streams, stream-ordered
allocation, zero-copy uploads, device pointers and pointer arithmetic,
bindless textures, cross-stream synchronisation, indirect launches driven by a
device-side count, graph capture and argument patching, and host callbacks.

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
   requires, so this one C struct is correct on every backend.
   --------------------------------------------------------------------------- */

typedef struct {
    md_gpu_uint4 dim;
    float        freq;
    uint32_t     _pad;
    uint64_t     vol;
} eval_args_t;

typedef struct {
    md_gpu_uint4 dim;
    float        threshold;
    uint32_t     _pad0;
    uint64_t     vol;
    uint64_t     count;
    uint64_t     indices;
    uint32_t     capacity;
    uint32_t     _pad1;
} compact_args_t;

typedef struct {
    md_gpu_uint4 dim;
    uint64_t     vol;
    uint64_t     count;
    uint64_t     indices;
    uint64_t     values;
} gather_args_t;

#define DEV(p) ((uint64_t)(uintptr_t)(p))

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
    uint32_t         count;
    float*           values;      /* host memory, filled by the callback */
    int              fired;
    md_gpu_ptr_t    staging;
    md_gpu_stream_t stream;
} readback_t;

static void on_readback_complete(void* user) {
    readback_t* rb = (readback_t*)user;
    rb->fired = 1;
    /* Everything the stream had issued before md_gpu_launch_host_fn has
       completed, so rb->count and rb->values are populated by now. */
    printf("  [callback] readback complete: %u voxels above threshold, "
           "first value %.4f\n", rb->count, rb->count ? rb->values[0] : 0.0f);
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
    /* One pool per memory kind -- a pool maps 1:1 onto the backend's memory
       object only if its kind is fixed at creation. */
    md_gpu_pool_desc_t pd = {0};
    pd.flags = MD_GPU_MEM_DEVICE;
    pd.label = "volume example";
    md_gpu_pool_t pool = md_gpu_pool_create(dev, &pd);
    CHECK(pool != NULL, "pool creation");

    /* -----------------------------------------------------------------
       Kernels
       ----------------------------------------------------------------- */
    md_gpu_kernel_desc_t kd = {0};
    kd.code = md_shader_gpu_volume_eval_field_start;
    kd.code_size = md_shader_gpu_volume_eval_field_byte_size;
    kd.label = "eval_field";
    md_gpu_kernel_t k_eval = md_gpu_kernel_create(dev, &kd);

    kd.code = md_shader_gpu_volume_compact_above_start;
    kd.code_size = md_shader_gpu_volume_compact_above_byte_size;
    kd.label = "compact_above";
    md_gpu_kernel_t k_compact = md_gpu_kernel_create(dev, &kd);

    kd.code = md_shader_gpu_volume_gather_start;
    kd.code_size = md_shader_gpu_volume_gather_byte_size;
    kd.label = "gather";
    md_gpu_kernel_t k_gather = md_gpu_kernel_create(dev, &kd);
    CHECK(k_eval && k_compact && k_gather, "kernel creation");

    md_gpu_kernel_info_t ki;
    md_gpu_kernel_info(k_eval, &ki);
    printf("eval_field threadgroup: %ux%ux%u\n",
           ki.group_size[0], ki.group_size[1], ki.group_size[2]);

    /* -----------------------------------------------------------------
       Resources. The texture handle goes straight into an argument struct;
       there is no binding and nothing to declare.
       ----------------------------------------------------------------- */
    md_gpu_tex_desc_t td = {0};
    td.width = DIM; td.height = DIM; td.depth = DIM;
    td.format = MD_GPU_FORMAT_R32_FLOAT;
    td.flags  = MD_GPU_TEX_STORAGE;
    td.label  = "field";
    md_gpu_tex_t vol = md_gpu_tex_create(dev, &td);
    CHECK(vol != 0, "texture creation");

    /* Stream-ordered allocation. Freeing these is legal at any point, even
       with work in flight, and the pool recycles them without a fence. */
    uint32_t* count   = md_gpu_malloc(pool, sizeof(uint32_t), compute);
    uint32_t* grid    = md_gpu_malloc(pool, 3 * sizeof(uint32_t), compute);
    uint32_t* indices = md_gpu_malloc(pool, VOXELS * sizeof(uint32_t), compute);
    float*    values  = md_gpu_malloc(pool, VOXELS * sizeof(float), compute);
    CHECK(count && grid && indices && values, "allocation");

    /* Device memory is a real pointer: arithmetic works and means what it
       means. `indices + 1024` is the 1025th element, no offset bookkeeping. */
    CHECK(md_gpu_ptr_base(indices + 1024) == (md_gpu_ptr_t)indices, "pointer arithmetic");

    /* A zero-copy upload, for data you would otherwise pack into a scratch
       buffer and memcpy. Here it seeds a lookup table nothing else uses,
       purely to show the shape. */
    float* seed = md_gpu_upload_begin(transfer, values, 16 * sizeof(float));
    CHECK(seed != NULL, "upload_begin");
    for (int i = 0; i < 16; ++i) seed[i] = (float)i;
    CHECK(md_gpu_upload_end(transfer), "upload_end");

    /* Join the two streams: compute waits for the transfer to land. */
    md_gpu_sync_t uploaded = md_gpu_stream_record(transfer);
    md_gpu_stream_wait(compute, uploaded);

    /* -----------------------------------------------------------------
       The pipeline. No barriers, no usage flags, no fences: everything
       issued into `compute` runs in order and sees the previous writes.
       ----------------------------------------------------------------- */
    const float threshold = 0.5f;

    eval_args_t ea = {0};
    ea.dim = (md_gpu_uint4){DIM, DIM, DIM, 0};
    ea.freq  = 4.0f;
    ea.vol   = vol;
    CHECK(md_gpu_launch(compute, k_eval, md_gpu_grid(DIM/4, DIM/4, DIM/4),
                         &ea, sizeof(ea)), "eval_field");

    CHECK(md_gpu_memset_async(count, 0, sizeof(uint32_t), compute), "memset count");

    compact_args_t ca = {0};
    ca.dim = (md_gpu_uint4){DIM, DIM, DIM, 0};
    ca.threshold = threshold;
    ca.vol       = vol;
    ca.count     = DEV(count);
    ca.indices   = DEV(indices);
    ca.capacity  = VOXELS;
    CHECK(md_gpu_launch(compute, k_compact, md_gpu_grid(DIM/4, DIM/4, DIM/4),
                         &ca, sizeof(ca)), "compact_above");

    /* Turn the device-side count into an indirect grid. The count is never
       read back to decide how much work to launch. */
    const uint32_t local[3] = {64, 1, 1};
    CHECK(md_gpu_make_grid(compute, grid, count, local), "make_grid");

    gather_args_t ga = {0};
    ga.dim = (md_gpu_uint4){DIM, DIM, DIM, 0};
    ga.vol     = vol;
    ga.count   = DEV(count);
    ga.indices = DEV(indices);
    ga.values  = DEV(values);
    CHECK(md_gpu_launch_indirect(compute, k_gather, grid,
                                  &ga, sizeof(ga)), "gather");

    /* -----------------------------------------------------------------
       Non-blocking readback. Nothing here waits.
       ----------------------------------------------------------------- */
    static float host_values[VOXELS];
    static readback_t rb;
    rb.values = host_values;
    rb.stream = compute;

    CHECK(md_gpu_memcpy_async(&rb.count, count, sizeof(uint32_t), compute), "readback count");
    CHECK(md_gpu_memcpy_async(host_values, values, sizeof(host_values), compute), "readback values");
    CHECK(md_gpu_launch_host_fn(compute, on_readback_complete, &rb), "host fn");

    /* In a real frame loop this is the only call you make; the callback fires
       on this thread, so touching OpenGL or ImGui from it is legal. */
    while (!rb.fired) {
        md_gpu_device_poll(dev);
    }

    /* Verify against the analytic answer. */
    uint32_t expect = 0;
    for (int z = 0; z < DIM; ++z)
    for (int y = 0; y < DIM; ++y)
    for (int x = 0; x < DIM; ++x) {
        float px = ((float)x + 0.5f) / DIM * 2.0f - 1.0f;
        float py = ((float)y + 0.5f) / DIM * 2.0f - 1.0f;
        float pz = ((float)z + 0.5f) / DIM * 2.0f - 1.0f;
        if (expf(-4.0f * (px*px + py*py + pz*pz)) > threshold) expect++;
    }
    printf("voxels above %.2f: gpu=%u cpu=%u %s\n",
           threshold, rb.count, expect, rb.count == expect ? "(match)" : "(MISMATCH)");
    CHECK(rb.count == expect, "count matches the analytic result");
    for (uint32_t i = 0; i < rb.count; ++i) {
        CHECK(host_values[i] > threshold, "every gathered value is above the threshold");
    }

    /* -----------------------------------------------------------------
       The same chain as a graph: record once, replay with new parameters.
       This is what a slider drag should cost.
       ----------------------------------------------------------------- */
    CHECK(md_gpu_capture_begin(compute, "volume pipeline"), "capture_begin");

    uint32_t eval_idx = md_gpu_capture_next_index(compute);
    md_gpu_launch(compute, k_eval, md_gpu_grid(DIM/4, DIM/4, DIM/4), &ea, sizeof(ea));
    md_gpu_memset_async(count, 0, sizeof(uint32_t), compute);
    uint32_t compact_idx = md_gpu_capture_next_index(compute);
    md_gpu_launch(compute, k_compact, md_gpu_grid(DIM/4, DIM/4, DIM/4), &ca, sizeof(ca));
    md_gpu_make_grid(compute, grid, count, local);
    md_gpu_launch_indirect(compute, k_gather, grid, &ga, sizeof(ga));

    md_gpu_graph_t graph = md_gpu_capture_end(compute);
    CHECK(graph != NULL, "capture_end");
    printf("graph captured: %u launches\n", md_gpu_graph_launch_count(graph));

    /* Patch the threshold in place and relaunch. The graph is not re-recorded
       and no command buffer is rebuilt -- this is a struct write. */
    const float thresholds[] = {0.25f, 0.75f};
    for (int t = 0; t < 2; ++t) {
        compact_args_t* patch = md_gpu_graph_args(graph, compact_idx);
        CHECK(patch != NULL, "graph_args");
        patch->threshold = thresholds[t];

        CHECK(md_gpu_graph_launch(graph, compute), "graph_launch");

        uint32_t n = 0;
        md_gpu_memcpy_async(&n, count, sizeof(n), compute);
        md_gpu_stream_sync(compute);

        uint32_t want = 0;
        for (int z = 0; z < DIM; ++z)
        for (int y = 0; y < DIM; ++y)
        for (int x = 0; x < DIM; ++x) {
            float px = ((float)x + 0.5f) / DIM * 2.0f - 1.0f;
            float py = ((float)y + 0.5f) / DIM * 2.0f - 1.0f;
            float pz = ((float)z + 0.5f) / DIM * 2.0f - 1.0f;
            if (expf(-4.0f * (px*px + py*py + pz*pz)) > thresholds[t]) want++;
        }
        printf("replay threshold %.2f: gpu=%u cpu=%u %s\n",
               thresholds[t], n, want, n == want ? "(match)" : "(MISMATCH)");
        CHECK(n == want, "graph replay honours the patched threshold");
    }
    (void)eval_idx;

    /* -----------------------------------------------------------------
       Teardown. Freeing is stream-ordered and never blocks; destroying the
       texture is legal even with work still in flight.
       ----------------------------------------------------------------- */
    md_gpu_pool_stats_t st;
    md_gpu_pool_stats(pool, &st);
    printf("pool: %llu in use, %llu reserved, peak %llu, %llu allocs (%llu from cache)\n",
           (unsigned long long)st.bytes_in_use, (unsigned long long)st.bytes_reserved,
           (unsigned long long)st.bytes_peak_in_use,
           (unsigned long long)st.alloc_count, (unsigned long long)st.reuse_count);

    md_gpu_graph_destroy(graph);

    /* One call instead of four md_gpu_free's. Stream-ordered, so it is safe
       even though the last graph launch may still be running. */
    md_gpu_pool_reset(pool, compute);
    md_gpu_tex_destroy(vol, compute);
    md_gpu_stream_sync(compute);
    md_gpu_device_poll(dev);

    md_gpu_pool_stats(pool, &st);
    printf("after reset: %llu in use, %llu still reserved for reuse\n",
           (unsigned long long)st.bytes_in_use, (unsigned long long)st.bytes_reserved);

    md_gpu_pool_destroy(pool);
    md_gpu_kernel_destroy(k_eval);
    md_gpu_kernel_destroy(k_compact);
    md_gpu_kernel_destroy(k_gather);
    md_gpu_device_destroy(dev);      /* waits for idle internally */

    printf("OK\n");
    return 0;
}
