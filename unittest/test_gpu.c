#include "utest.h"

#include <core/md_gpu.h>
#include <core/md_allocator.h>
#include <core/md_os.h>

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#if MD_ENABLE_GPU

#include "gpu_test_shaders.inl"

/* =========================================================================
   Fixtures
   ========================================================================= */

/* A counting allocator, used to prove md_gpu routes every host-side
   allocation through md_gpu_device_desc_t::alloc and releases all of it. */
typedef struct {
    size_t live_bytes;
    size_t total_bytes;
    size_t alloc_count;
} gpu_alloc_stats_t;

static void* gpu_test_realloc(md_allocator_o* inst, void* ptr, size_t old_size, size_t new_size, const char* file, size_t line) {
    (void)file; (void)line;
    gpu_alloc_stats_t* stats = (gpu_alloc_stats_t*)inst;
    stats->live_bytes -= old_size;
    if (new_size == 0) {
        free(ptr);
        return NULL;
    }
    void* mem = realloc(ptr, new_size);
    if (!mem) return NULL;
    stats->live_bytes  += new_size;
    stats->total_bytes += new_size;
    stats->alloc_count += 1;
    return mem;
}

typedef struct {
    md_gpu_device_t dev;
    md_gpu_stream_t compute;
    md_gpu_stream_t transfer;
    md_gpu_pool_t   pool;       /* device-local  */
    md_gpu_pool_t   read_pool;  /* host-readable */
    md_gpu_kernel_t k_fill;
    md_gpu_kernel_t k_scale;
    md_gpu_kernel_t k_sum;
    md_gpu_kernel_t k_tex_write;
    md_gpu_kernel_t k_tex_read;
    md_gpu_kernel_t k_bump;
    md_gpu_kernel_t k_layout;
    md_gpu_kernel_t k_tex_probe;
} gpu_fixture_t;

/* group_size is NOT optional. Vulkan recovers it from the SPIR-V when left
   zero, but Metal cannot read it out of a metallib, so md_gpu_kernel_create
   falls back to {1,1,1} there -- which turns md_gpu_grid_1d(N, 64) into N/64
   threads instead of N and makes every data-carrying test fail for a reason
   that has nothing to do with what it is testing. Declare it always; it must
   match [numthreads(...)] in gpu_test.slang. */
#define GPU_KERNEL(fix, field, sym, name, gx, gy, gz)                             \
    do {                                                                          \
        md_gpu_kernel_desc_t kd = {0};                                            \
        kd.code          = sym##_start;                                           \
        kd.code_size     = sym##_byte_size;                                       \
        kd.label         = name;                                                  \
        kd.group_size[0] = (gx);                                                  \
        kd.group_size[1] = (gy);                                                  \
        kd.group_size[2] = (gz);                                                  \
        (fix)->field = md_gpu_kernel_create((fix)->dev, &kd);                     \
    } while (0)

static bool gpu_open(gpu_fixture_t* f) {
    memset(f, 0, sizeof(*f));
    md_gpu_device_desc_t dd = {0};
    dd.enable_validation = true;
    dd.label             = "md_gpu unittest";
    f->dev = md_gpu_device_create(&dd);
    if (!f->dev) return false;

    f->compute  = md_gpu_stream_default(f->dev, MD_GPU_STREAM_COMPUTE);
    f->transfer = md_gpu_stream_default(f->dev, MD_GPU_STREAM_TRANSFER);
    md_gpu_pool_desc_t pd = {0};
    pd.flags = MD_GPU_MEM_DEVICE;   pd.label = "test device";   f->pool      = md_gpu_pool_create(f->dev, &pd);
    pd.flags = MD_GPU_MEM_HOST_READ; pd.label = "test readback"; f->read_pool = md_gpu_pool_create(f->dev, &pd);
    if (!f->pool || !f->read_pool) return false;

    GPU_KERNEL(f, k_fill,      md_shader_gpu_test_fill,         "fill",         64, 1, 1);
    GPU_KERNEL(f, k_scale,     md_shader_gpu_test_scale_add,    "scale_add",    64, 1, 1);
    GPU_KERNEL(f, k_sum,       md_shader_gpu_test_sum_reduce,   "sum_reduce",   64, 1, 1);
    GPU_KERNEL(f, k_tex_write, md_shader_gpu_test_tex_write,    "tex_write",     4, 4, 4);
    GPU_KERNEL(f, k_tex_read,  md_shader_gpu_test_tex_read,     "tex_read",      4, 4, 4);
    GPU_KERNEL(f, k_bump,      md_shader_gpu_test_bump,         "bump",         64, 1, 1);
    GPU_KERNEL(f, k_layout,    md_shader_gpu_test_layout_probe, "layout_probe",  1, 1, 1);
    GPU_KERNEL(f, k_tex_probe, md_shader_gpu_test_tex_probe,    "tex_probe",     1, 1, 1);
    return f->k_fill && f->k_scale && f->k_sum && f->k_tex_write &&
           f->k_tex_read && f->k_bump && f->k_layout && f->k_tex_probe;
}

/* Why the fixture could not be opened -- no Vulkan loader, no driver (ICD), no
   compute queue, a kernel that failed to build. Without this a CI log shows
   only "skipped" and gives no way to tell a missing GPU from a broken build. */
static const char* gpu_no_device_reason(void) {
    const char* err = md_gpu_last_error();
    return (err && err[0]) ? err : "No GPU device available";
}

static void gpu_close(gpu_fixture_t* f) {
    if (f->k_fill)      md_gpu_kernel_destroy(f->k_fill);
    if (f->k_scale)     md_gpu_kernel_destroy(f->k_scale);
    if (f->k_sum)       md_gpu_kernel_destroy(f->k_sum);
    if (f->k_tex_write) md_gpu_kernel_destroy(f->k_tex_write);
    if (f->k_tex_read)  md_gpu_kernel_destroy(f->k_tex_read);
    if (f->k_bump)      md_gpu_kernel_destroy(f->k_bump);
    if (f->k_layout)    md_gpu_kernel_destroy(f->k_layout);
    if (f->k_tex_probe) md_gpu_kernel_destroy(f->k_tex_probe);
    if (f->pool)        md_gpu_pool_destroy(f->pool);
    if (f->read_pool)   md_gpu_pool_destroy(f->read_pool);
    md_gpu_device_destroy(f->dev);
}

/* Argument structs, mirroring unittest/shaders/gpu_test.slang. */
typedef struct { uint32_t n, base, pad0, pad1; uint64_t dst; }                 fill_args_t;
typedef struct { uint32_t n, mul, add, pad; uint64_t src, dst; }               scale_args_t;
typedef struct { uint32_t n, pad0, pad1, pad2; uint64_t src, out_val; }        sum_args_t;
typedef struct { uint32_t dim[3]; float scale; uint64_t tex; }                 tex_args_t;
typedef struct { uint32_t dim[3], pad; uint64_t tex, dst; }                    tex_read_args_t;
typedef struct { uint32_t n, delta, pad0, pad1; uint64_t dst; }                bump_args_t;

/* The two probe structs are the point of the exercise, so they are mirrored
   with md_gpu.h's vector types rather than raw arrays -- that is the machinery
   under test. */
typedef struct {
    md_gpu_float4   v4;
    uint32_t        dim_x;
    uint32_t        dim_y;
    md_gpu_uint2    pair;
    uint32_t        dim_z;
    float           scale;
    uint64_t        dst;
} layout_probe_args_t;

typedef struct {
    uint32_t     n, pad0, pad1, pad2;
    uint64_t     dst;
    md_gpu_tex_t tex;
    uint32_t     marker, pad3;
} tex_probe_args_t;

#define DEV_ADDR(p) ((uint64_t)(uintptr_t)(p))

/* =========================================================================
   Device and streams
   ========================================================================= */

UTEST(gpu, device_create_destroy) {
    md_gpu_device_t dev = md_gpu_device_create(NULL);
    if (!dev) UTEST_SKIP(gpu_no_device_reason());
    md_gpu_device_info_t info;
    ASSERT_TRUE(md_gpu_device_info(dev, &info));
    ASSERT_TRUE(info.name[0] != '\0');
    ASSERT_GT(info.max_threads_per_group, 0u);
    md_gpu_device_destroy(dev);
}

UTEST(gpu, device_routes_all_host_allocations) {
    gpu_alloc_stats_t stats = {0};
    md_allocator_i alloc = { (md_allocator_o*)&stats, gpu_test_realloc };

    md_gpu_device_desc_t dd = {0};
    dd.alloc = &alloc;
    md_gpu_device_t dev = md_gpu_device_create(&dd);
    if (!dev) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_stream_t s = md_gpu_stream_default(dev, MD_GPU_STREAM_COMPUTE);
    md_gpu_pool_desc_t pd = {0};
    pd.flags = MD_GPU_MEM_DEVICE;
    md_gpu_pool_t pool = md_gpu_pool_create(dev, &pd);
    ASSERT_TRUE(pool != NULL);
    md_gpu_ptr_t p = md_gpu_malloc(pool, 4096, s);
    ASSERT_TRUE(p != NULL);
    md_gpu_free(p, s);
    md_gpu_pool_destroy(pool);

    ASSERT_GT(stats.alloc_count, 0u);
    md_gpu_device_destroy(dev);
    EXPECT_EQ(0u, (unsigned)stats.live_bytes);
}

UTEST(gpu, streams_create_and_sync) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_stream_t a = md_gpu_stream_create(f.dev, MD_GPU_STREAM_COMPUTE, "a");
    ASSERT_TRUE(a != NULL);
    ASSERT_TRUE(md_gpu_stream_device(a) == f.dev);

    /* A stream with nothing issued reports a none sync and syncs instantly. */
    md_gpu_sync_t none = md_gpu_stream_record(a);
    EXPECT_FALSE(md_gpu_sync_is_valid(none));
    EXPECT_TRUE(md_gpu_sync_is_complete(none));
    md_gpu_stream_sync(a);

    md_gpu_stream_destroy(a);
    gpu_close(&f);
}

/* =========================================================================
   Memory
   ========================================================================= */

UTEST(gpu, memcpy_roundtrip) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 1024 };
    uint32_t src[N], dst[N];
    for (int i = 0; i < N; ++i) src[i] = (uint32_t)i * 2654435761u;
    memset(dst, 0, sizeof(dst));

    uint32_t* d = (uint32_t*)md_gpu_malloc(f.pool, sizeof(src), f.compute);
    ASSERT_TRUE(d != NULL);
    ASSERT_TRUE(md_gpu_memcpy_async(d, src, sizeof(src), f.compute));
    ASSERT_TRUE(md_gpu_memcpy_async(dst, d, sizeof(dst), f.compute));
    md_gpu_stream_sync(f.compute);

    for (int i = 0; i < N; ++i) EXPECT_EQ(src[i], dst[i]);

    md_gpu_free(d, f.compute);
    gpu_close(&f);
}

UTEST(gpu, memcpy_device_to_device) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 256 };
    uint32_t src[N], dst[N];
    for (int i = 0; i < N; ++i) src[i] = (uint32_t)(i + 7) * 11u;
    memset(dst, 0, sizeof(dst));

    uint32_t* a = (uint32_t*)md_gpu_malloc(f.pool, sizeof(src), f.compute);
    uint32_t* b = (uint32_t*)md_gpu_malloc(f.pool, sizeof(src), f.compute);
    ASSERT_TRUE(a && b);

    ASSERT_TRUE(md_gpu_memcpy_async(a, src, sizeof(src), f.compute));
    ASSERT_TRUE(md_gpu_memcpy_async(b, a, sizeof(src), f.compute));
    ASSERT_TRUE(md_gpu_memcpy_async(dst, b, sizeof(dst), f.compute));
    md_gpu_stream_sync(f.compute);

    for (int i = 0; i < N; ++i) EXPECT_EQ(src[i], dst[i]);

    md_gpu_free(a, f.compute);
    md_gpu_free(b, f.compute);
    gpu_close(&f);
}

UTEST(gpu, memset_aligned_and_unaligned) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { BYTES = 256 };
    uint8_t host[BYTES];
    uint8_t* d = (uint8_t*)md_gpu_malloc(f.pool, BYTES, f.compute);
    ASSERT_TRUE(d != NULL);

    /* Whole buffer, 4-byte aligned. */
    ASSERT_TRUE(md_gpu_memset_async(d, 0xAB, BYTES, f.compute));
    memset(host, 0, BYTES);
    ASSERT_TRUE(md_gpu_memcpy_async(host, d, BYTES, f.compute));
    md_gpu_stream_sync(f.compute);
    for (int i = 0; i < BYTES; ++i) EXPECT_EQ(0xAB, host[i]);

    /* Unaligned offset and length, exercising the staged head/tail path. */
    ASSERT_TRUE(md_gpu_memset_async(d + 3, 0x5C, 10, f.compute));
    memset(host, 0, BYTES);
    ASSERT_TRUE(md_gpu_memcpy_async(host, d, BYTES, f.compute));
    md_gpu_stream_sync(f.compute);
    for (int i = 0; i < BYTES; ++i) {
        uint8_t expect = (i >= 3 && i < 13) ? 0x5C : 0xAB;
        EXPECT_EQ(expect, host[i]);
    }

    md_gpu_free(d, f.compute);
    gpu_close(&f);
}

UTEST(gpu, pointer_arithmetic_subranges) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 128 };
    uint32_t* base = (uint32_t*)md_gpu_malloc(f.pool, N * 2 * sizeof(uint32_t), f.compute);
    ASSERT_TRUE(base != NULL);

    /* The second half is addressed by plain pointer arithmetic. */
    uint32_t* upper = base + N;
    ASSERT_TRUE(md_gpu_ptr_base(upper) == (md_gpu_ptr_t)base);
    ASSERT_EQ((size_t)(N * sizeof(uint32_t)), md_gpu_ptr_size(upper));

    fill_args_t lo = {0}; lo.n = N; lo.base = 1000; lo.dst = DEV_ADDR(base);
    fill_args_t hi = {0}; hi.n = N; hi.base = 5000; hi.dst = DEV_ADDR(upper);
    ASSERT_TRUE(md_gpu_launch(f.compute, f.k_fill, md_gpu_grid_1d(N, 64), &lo, sizeof(lo)));
    ASSERT_TRUE(md_gpu_launch(f.compute, f.k_fill, md_gpu_grid_1d(N, 64), &hi, sizeof(hi)));

    uint32_t host[N * 2];
    ASSERT_TRUE(md_gpu_memcpy_async(host, base, sizeof(host), f.compute));
    md_gpu_stream_sync(f.compute);

    for (int i = 0; i < N; ++i) {
        EXPECT_EQ((uint32_t)(1000 + i), host[i]);
        EXPECT_EQ((uint32_t)(5000 + i), host[N + i]);
    }

    md_gpu_free(base, f.compute);
    gpu_close(&f);
}

UTEST(gpu, upload_begin_end_zero_copy) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 512 };
    uint32_t* d = (uint32_t*)md_gpu_malloc(f.pool, N * sizeof(uint32_t), f.compute);
    ASSERT_TRUE(d != NULL);

    uint32_t* p = (uint32_t*)md_gpu_upload_begin(f.compute, d, N * sizeof(uint32_t));
    ASSERT_TRUE(p != NULL);
    for (int i = 0; i < N; ++i) p[i] = (uint32_t)(i * 3 + 1);
    ASSERT_TRUE(md_gpu_upload_end(f.compute));

    uint32_t host[N];
    ASSERT_TRUE(md_gpu_memcpy_async(host, d, sizeof(host), f.compute));
    md_gpu_stream_sync(f.compute);
    for (int i = 0; i < N; ++i) EXPECT_EQ((uint32_t)(i * 3 + 1), host[i]);

    md_gpu_free(d, f.compute);
    gpu_close(&f);
}

UTEST(gpu, pool_reuses_freed_blocks) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_pool_t pool = md_gpu_pool_create(f.dev, &(md_gpu_pool_desc_t){ .flags = MD_GPU_MEM_DEVICE, .release_threshold = 16 * 1024 * 1024, .label = "reuse" });
    ASSERT_TRUE(pool != NULL);

    md_gpu_pool_stats_t st = {0};
    md_gpu_ptr_t a = md_gpu_malloc(pool, 4096, f.compute);
    ASSERT_TRUE(a != NULL);
    md_gpu_pool_stats(pool, &st);
    EXPECT_EQ(4096u, (unsigned)st.bytes_in_use);
    EXPECT_EQ(1u, st.blocks_in_use);
    EXPECT_EQ(1u, (unsigned)st.alloc_count);
    EXPECT_EQ(0u, (unsigned)st.reuse_count);
    uint64_t reserved = st.bytes_reserved;

    md_gpu_free(a, f.compute);
    md_gpu_pool_stats(pool, &st);
    EXPECT_EQ(0u, (unsigned)st.bytes_in_use);
    EXPECT_EQ(1u, st.blocks_cached);

    /* Freed on this stream, so reuse is immediate and hits the same block. */
    md_gpu_ptr_t b = md_gpu_malloc(pool, 4096, f.compute);
    EXPECT_TRUE(b == a);

    md_gpu_pool_stats(pool, &st);
    EXPECT_EQ(reserved, st.bytes_reserved);   /* no new memory was committed */
    EXPECT_EQ(1u, (unsigned)st.reuse_count);
    EXPECT_EQ(4096u, (unsigned)st.bytes_peak_in_use);

    md_gpu_free(b, f.compute);
    md_gpu_pool_destroy(pool);
    gpu_close(&f);
}

/* The CPU arena-reset pattern: drop every allocation in one call, keep the
   memory, and let the next round be served entirely from cache. */
UTEST(gpu, pool_reset_recycles_everything) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_pool_t pool = md_gpu_pool_create(f.dev, &(md_gpu_pool_desc_t){ .flags = MD_GPU_MEM_DEVICE, .release_threshold = 16 * 1024 * 1024, .label = "reset" });
    ASSERT_TRUE(pool != NULL);

    enum { N = 8 };
    md_gpu_ptr_t first[N];
    for (int i = 0; i < N; ++i) {
        first[i] = md_gpu_malloc(pool, 4096, f.compute);
        ASSERT_TRUE(first[i] != NULL);
    }

    md_gpu_pool_stats_t st = {0};
    md_gpu_pool_stats(pool, &st);
    EXPECT_EQ((unsigned)(N * 4096), (unsigned)st.bytes_in_use);
    EXPECT_EQ((uint32_t)N, st.blocks_in_use);
    uint64_t reserved_before = st.bytes_reserved;

    /* One call frees them all, without giving the memory back. */
    md_gpu_pool_reset(pool, f.compute);

    md_gpu_pool_stats(pool, &st);
    EXPECT_EQ(0u, (unsigned)st.bytes_in_use);
    EXPECT_EQ(0u, st.blocks_in_use);
    EXPECT_EQ((uint32_t)N, st.blocks_cached);
    EXPECT_EQ(reserved_before, st.bytes_reserved);      /* memory retained */

    /* The next round is served entirely from cache and commits nothing new. */
    uint64_t reuse_before = st.reuse_count;
    md_gpu_ptr_t second[N];
    for (int i = 0; i < N; ++i) {
        second[i] = md_gpu_malloc(pool, 4096, f.compute);
        ASSERT_TRUE(second[i] != NULL);
    }
    md_gpu_pool_stats(pool, &st);
    EXPECT_EQ(reserved_before, st.bytes_reserved);
    EXPECT_EQ(reuse_before + N, st.reuse_count);

    /* Reset again, then hand the memory back for real. */
    md_gpu_pool_reset(pool, f.compute);
    md_gpu_stream_sync(f.compute);
    md_gpu_pool_trim(pool, 0);
    md_gpu_pool_stats(pool, &st);
    EXPECT_EQ(0u, (unsigned)st.bytes_reserved);

    md_gpu_pool_destroy(pool);
    gpu_close(&f);
}

/* Reset while the GPU is still reading the memory must not disturb it: the
   blocks only become reusable at that point in the stream. */
UTEST(gpu, pool_reset_is_stream_ordered) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_pool_t pool = md_gpu_pool_create(f.dev, &(md_gpu_pool_desc_t){ .flags = MD_GPU_MEM_DEVICE, .release_threshold = 16 * 1024 * 1024, .label = "reset-order" });
    ASSERT_TRUE(pool != NULL);

    enum { N = 1024 };
    uint32_t* d = (uint32_t*)md_gpu_malloc(pool, N * sizeof(uint32_t), f.compute);
    ASSERT_TRUE(d != NULL);

    fill_args_t a = {0};
    a.n = N; a.base = 5; a.dst = DEV_ADDR(d);
    ASSERT_TRUE(md_gpu_launch(f.compute, f.k_fill, md_gpu_grid_1d(N, 64), &a, sizeof(a)));

    static uint32_t host[N];
    ASSERT_TRUE(md_gpu_memcpy_async(host, d, sizeof(host), f.compute));

    /* Reset with the launch and the readback still in flight. */
    md_gpu_pool_reset(pool, f.compute);

    md_gpu_stream_sync(f.compute);
    for (int i = 0; i < N; ++i) EXPECT_EQ((uint32_t)(5 + i), host[i]);

    md_gpu_pool_destroy(pool);
    gpu_close(&f);
}

UTEST(gpu, pool_trim_releases) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_pool_t pool = md_gpu_pool_create(f.dev, &(md_gpu_pool_desc_t){ .flags = MD_GPU_MEM_DEVICE, .release_threshold = 64 * 1024 * 1024, .label = "trim" });
    ASSERT_TRUE(pool != NULL);

    md_gpu_ptr_t p = md_gpu_malloc(pool, 1024 * 1024, f.compute);
    ASSERT_TRUE(p != NULL);
    md_gpu_free(p, f.compute);
    md_gpu_stream_sync(f.compute);

    md_gpu_pool_stats_t st = {0};
    md_gpu_pool_stats(pool, &st);
    EXPECT_GT(st.bytes_reserved, 0u);

    md_gpu_pool_trim(pool, 0);
    md_gpu_pool_stats(pool, &st);
    EXPECT_EQ(0u, (unsigned)st.bytes_reserved);

    md_gpu_pool_destroy(pool);
    gpu_close(&f);
}

/* =========================================================================
   Launches and program order
   ========================================================================= */

UTEST(gpu, launch_writes_buffer) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 1000 };
    uint32_t* d = (uint32_t*)md_gpu_malloc(f.pool, N * sizeof(uint32_t), f.compute);
    ASSERT_TRUE(d != NULL);

    fill_args_t a = {0};
    a.n = N; a.base = 100; a.dst = DEV_ADDR(d);
    ASSERT_TRUE(md_gpu_launch(f.compute, f.k_fill, md_gpu_grid_1d(N, 64), &a, sizeof(a)));

    uint32_t host[N];
    ASSERT_TRUE(md_gpu_memcpy_async(host, d, sizeof(host), f.compute));
    md_gpu_stream_sync(f.compute);

    for (int i = 0; i < N; ++i) EXPECT_EQ((uint32_t)(100 + i), host[i]);

    md_gpu_free(d, f.compute);
    gpu_close(&f);
}

/* The central guarantee: consecutive launches in one stream see each other's
   writes with no barrier, no usage declaration and no fence from the caller. */
UTEST(gpu, program_order_chain) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 4096, STEPS = 40 };
    uint32_t* a = (uint32_t*)md_gpu_malloc(f.pool, N * sizeof(uint32_t), f.compute);
    uint32_t* b = (uint32_t*)md_gpu_malloc(f.pool, N * sizeof(uint32_t), f.compute);
    ASSERT_TRUE(a && b);

    fill_args_t fa = {0};
    fa.n = N; fa.base = 0; fa.dst = DEV_ADDR(a);
    ASSERT_TRUE(md_gpu_launch(f.compute, f.k_fill, md_gpu_grid_1d(N, 64), &fa, sizeof(fa)));

    /* Ping-pong a -> b -> a -> ... far past the old MD_GPU_MAX_DISPATCHES_PER_CMD
       limit of 32, proving there is no per-command-buffer dispatch wall. */
    uint32_t* src = a;
    uint32_t* dst = b;
    for (int step = 0; step < STEPS; ++step) {
        scale_args_t sa = {0};
        sa.n = N; sa.mul = 1; sa.add = 1;
        sa.src = DEV_ADDR(src);
        sa.dst = DEV_ADDR(dst);
        ASSERT_TRUE(md_gpu_launch(f.compute, f.k_scale, md_gpu_grid_1d(N, 64), &sa, sizeof(sa)));
        uint32_t* tmp = src; src = dst; dst = tmp;
    }

    uint32_t host[N];
    ASSERT_TRUE(md_gpu_memcpy_async(host, src, sizeof(host), f.compute));
    md_gpu_stream_sync(f.compute);

    for (int i = 0; i < N; ++i) EXPECT_EQ((uint32_t)(i + STEPS), host[i]);

    md_gpu_free(a, f.compute);
    md_gpu_free(b, f.compute);
    gpu_close(&f);
}

/* The same chain, but split across several submissions: the dependency must
   survive command-buffer and submission boundaries. */
UTEST(gpu, program_order_across_submissions) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 1024, STEPS = 8 };
    uint32_t* a = (uint32_t*)md_gpu_malloc(f.pool, N * sizeof(uint32_t), f.compute);
    uint32_t* b = (uint32_t*)md_gpu_malloc(f.pool, N * sizeof(uint32_t), f.compute);
    ASSERT_TRUE(a && b);

    fill_args_t fa = {0};
    fa.n = N; fa.base = 0; fa.dst = DEV_ADDR(a);
    ASSERT_TRUE(md_gpu_launch(f.compute, f.k_fill, md_gpu_grid_1d(N, 64), &fa, sizeof(fa)));

    uint32_t* src = a;
    uint32_t* dst = b;
    for (int step = 0; step < STEPS; ++step) {
        scale_args_t sa = {0};
        sa.n = N; sa.mul = 2; sa.add = 0;
        sa.src = DEV_ADDR(src);
        sa.dst = DEV_ADDR(dst);
        ASSERT_TRUE(md_gpu_launch(f.compute, f.k_scale, md_gpu_grid_1d(N, 64), &sa, sizeof(sa)));
        md_gpu_stream_flush(f.compute);          /* force a submission boundary */
        uint32_t* tmp = src; src = dst; dst = tmp;
    }

    uint32_t host[N];
    ASSERT_TRUE(md_gpu_memcpy_async(host, src, sizeof(host), f.compute));
    md_gpu_stream_sync(f.compute);

    for (int i = 0; i < N; ++i) EXPECT_EQ(((uint32_t)i << STEPS), host[i]);

    md_gpu_free(a, f.compute);
    md_gpu_free(b, f.compute);
    gpu_close(&f);
}

UTEST(gpu, cross_stream_dependency) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 2048 };
    uint32_t src[N];
    for (int i = 0; i < N; ++i) src[i] = (uint32_t)i;

    uint32_t* d   = (uint32_t*)md_gpu_malloc(f.pool, sizeof(src), f.transfer);
    uint32_t* out = (uint32_t*)md_gpu_malloc(f.pool, sizeof(src), f.transfer);
    ASSERT_TRUE(d && out);

    /* Upload on the transfer stream. */
    ASSERT_TRUE(md_gpu_memcpy_async(d, src, sizeof(src), f.transfer));
    md_gpu_sync_t uploaded = md_gpu_stream_record(f.transfer);
    ASSERT_TRUE(md_gpu_sync_is_valid(uploaded));

    /* Compute stream waits for it, then consumes the data. */
    md_gpu_stream_wait(f.compute, uploaded);
    scale_args_t sa = {0};
    sa.n = N; sa.mul = 3; sa.add = 5;
    sa.src = DEV_ADDR(d);
    sa.dst = DEV_ADDR(out);
    ASSERT_TRUE(md_gpu_launch(f.compute, f.k_scale, md_gpu_grid_1d(N, 64), &sa, sizeof(sa)));

    uint32_t host[N];
    ASSERT_TRUE(md_gpu_memcpy_async(host, out, sizeof(host), f.compute));
    md_gpu_stream_sync(f.compute);

    for (int i = 0; i < N; ++i) EXPECT_EQ((uint32_t)(i * 3 + 5), host[i]);

    md_gpu_free(d, f.transfer);
    md_gpu_free(out, f.transfer);
    gpu_close(&f);
}

/* Waiting on a sync from your own stream is a no-op, and a none sync is too. */
UTEST(gpu, wait_on_none_and_self_is_noop) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 64 };
    uint32_t* d = (uint32_t*)md_gpu_malloc(f.pool, N * sizeof(uint32_t), f.compute);
    ASSERT_TRUE(d != NULL);

    md_gpu_stream_wait(f.compute, md_gpu_sync_none());

    fill_args_t a = {0};
    a.n = N; a.base = 42; a.dst = DEV_ADDR(d);
    ASSERT_TRUE(md_gpu_launch(f.compute, f.k_fill, md_gpu_grid_1d(N, 64), &a, sizeof(a)));

    md_gpu_sync_t self = md_gpu_stream_record(f.compute);
    md_gpu_stream_wait(f.compute, self);          /* must not deadlock */

    uint32_t host[N];
    ASSERT_TRUE(md_gpu_memcpy_async(host, d, sizeof(host), f.compute));
    md_gpu_stream_sync(f.compute);
    for (int i = 0; i < N; ++i) EXPECT_EQ((uint32_t)(42 + i), host[i]);

    md_gpu_free(d, f.compute);
    gpu_close(&f);
}

/* =========================================================================
   Textures
   ========================================================================= */

UTEST(gpu, texture_kernel_write_then_host_read) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { D = 8, VOXELS = D * D * D };
    md_gpu_tex_desc_t td = {0};
    td.width = D; td.height = D; td.depth = D;
    td.format = MD_GPU_FORMAT_R32_FLOAT;
    td.flags  = MD_GPU_TEX_STORAGE;
    td.label  = "volume";
    md_gpu_tex_t tex = md_gpu_tex_create(f.dev, &td);
    ASSERT_TRUE(tex != 0);

    md_gpu_tex_desc_t back = {0};
    ASSERT_TRUE(md_gpu_tex_desc(tex, &back));
    EXPECT_EQ((uint32_t)D, back.width);
    EXPECT_EQ((uint32_t)D, back.depth);

    tex_args_t ta = {0};
    ta.dim[0] = D; ta.dim[1] = D; ta.dim[2] = D;
    ta.tex = tex;
    ta.scale = 2.0f;
    ASSERT_TRUE(md_gpu_launch(f.compute, f.k_tex_write, md_gpu_grid(D/4, D/4, D/4), &ta, sizeof(ta)));

    float host[VOXELS];
    memset(host, 0, sizeof(host));
    ASSERT_TRUE(md_gpu_memcpy_from_tex_async(host, tex, NULL, sizeof(host), f.compute));
    md_gpu_stream_sync(f.compute);

    for (int i = 0; i < VOXELS; ++i) EXPECT_EQ((float)i * 2.0f, host[i]);

    md_gpu_tex_destroy(tex, f.compute);
    md_gpu_device_poll(f.dev);
    gpu_close(&f);
}

UTEST(gpu, texture_host_write_then_kernel_read) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { D = 8, VOXELS = D * D * D };
    md_gpu_tex_desc_t td = {0};
    td.width = D; td.height = D; td.depth = D;
    td.format = MD_GPU_FORMAT_R32_FLOAT;
    td.flags  = MD_GPU_TEX_STORAGE;
    md_gpu_tex_t tex = md_gpu_tex_create(f.dev, &td);
    ASSERT_TRUE(tex != 0);

    float src[VOXELS];
    for (int i = 0; i < VOXELS; ++i) src[i] = (float)(VOXELS - i);
    ASSERT_TRUE(md_gpu_memcpy_to_tex_async(tex, NULL, src, sizeof(src), f.compute));

    float* d = (float*)md_gpu_malloc(f.pool, sizeof(src), f.compute);
    ASSERT_TRUE(d != NULL);

    tex_read_args_t ra = {0};
    ra.dim[0] = D; ra.dim[1] = D; ra.dim[2] = D;
    ra.tex = tex;
    ra.dst = DEV_ADDR(d);
    ASSERT_TRUE(md_gpu_launch(f.compute, f.k_tex_read, md_gpu_grid(D/4, D/4, D/4), &ra, sizeof(ra)));

    float host[VOXELS];
    memset(host, 0, sizeof(host));
    ASSERT_TRUE(md_gpu_memcpy_async(host, d, sizeof(host), f.compute));
    md_gpu_stream_sync(f.compute);

    for (int i = 0; i < VOXELS; ++i) EXPECT_EQ(src[i], host[i]);

    md_gpu_free(d, f.compute);
    md_gpu_tex_destroy(tex, f.compute);
    md_gpu_device_poll(f.dev);
    gpu_close(&f);
}

/* Destroying a texture while work that used it is still in flight must be
   safe, and the slot must come back after a poll. */
UTEST(gpu, texture_destroy_while_in_flight) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { D = 16 };
    md_gpu_tex_desc_t td = {0};
    td.width = D; td.height = D; td.depth = D;
    td.format = MD_GPU_FORMAT_R32_FLOAT;
    td.flags  = MD_GPU_TEX_STORAGE;

    for (int iter = 0; iter < 8; ++iter) {
        md_gpu_tex_t tex = md_gpu_tex_create(f.dev, &td);
        ASSERT_TRUE(tex != 0);
        tex_args_t ta = {0};
        ta.dim[0] = D; ta.dim[1] = D; ta.dim[2] = D;
        ta.tex = tex;
        ta.scale = 1.0f;
        ASSERT_TRUE(md_gpu_launch(f.compute, f.k_tex_write, md_gpu_grid(D/4, D/4, D/4), &ta, sizeof(ta)));
        md_gpu_stream_flush(f.compute);
        md_gpu_tex_destroy(tex, f.compute);       /* still in flight — legal */
        md_gpu_device_poll(f.dev);
    }

    md_gpu_stream_sync(f.compute);
    md_gpu_device_poll(f.dev);
    gpu_close(&f);
}

/* =========================================================================
   Indirect launches
   ========================================================================= */

UTEST(gpu, make_grid_and_indirect_launch) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 500 };

    uint32_t* data  = (uint32_t*)md_gpu_malloc(f.pool, N * sizeof(uint32_t), f.compute);
    uint32_t* ones  = (uint32_t*)md_gpu_malloc(f.pool, N * sizeof(uint32_t), f.compute);
    uint32_t* count = (uint32_t*)md_gpu_malloc(f.pool, sizeof(uint32_t), f.compute);
    uint32_t* grid  = (uint32_t*)md_gpu_malloc(f.pool, 3 * sizeof(uint32_t), f.compute);
    ASSERT_TRUE(data && ones && count && grid);

    /* data[i] = i */
    fill_args_t fa = {0};
    fa.n = N; fa.base = 0; fa.dst = DEV_ADDR(data);
    ASSERT_TRUE(md_gpu_launch(f.compute, f.k_fill, md_gpu_grid_1d(N, 64), &fa, sizeof(fa)));

    /* ones[i] = 1, so the reduction produces exactly N on the device. */
    scale_args_t za = {0};
    za.n = N; za.mul = 0; za.add = 1;
    za.src = DEV_ADDR(data);
    za.dst = DEV_ADDR(ones);
    ASSERT_TRUE(md_gpu_launch(f.compute, f.k_scale, md_gpu_grid_1d(N, 64), &za, sizeof(za)));

    ASSERT_TRUE(md_gpu_memset_async(count, 0, sizeof(uint32_t), f.compute));

    sum_args_t sa = {0};
    sa.n = N; sa.src = DEV_ADDR(ones); sa.out_val = DEV_ADDR(count);
    ASSERT_TRUE(md_gpu_launch(f.compute, f.k_sum, md_gpu_grid_1d(N, 64), &sa, sizeof(sa)));

    /* The count never touches the CPU: it drives the indirect grid directly. */
    const uint32_t local[3] = {64, 1, 1};
    ASSERT_TRUE(md_gpu_make_grid(f.compute, grid, count, local));

    bump_args_t ba = {0};
    ba.n = N; ba.delta = 1000; ba.dst = DEV_ADDR(data);
    ASSERT_TRUE(md_gpu_launch_indirect(f.compute, f.k_bump, grid, &ba, sizeof(ba)));

    uint32_t host_count = 0, host_grid[3] = {0}, host[N];
    ASSERT_TRUE(md_gpu_memcpy_async(&host_count, count, sizeof(host_count), f.compute));
    ASSERT_TRUE(md_gpu_memcpy_async(host_grid, grid, sizeof(host_grid), f.compute));
    ASSERT_TRUE(md_gpu_memcpy_async(host, data, sizeof(host), f.compute));
    md_gpu_stream_sync(f.compute);

    EXPECT_EQ((uint32_t)N, host_count);
    EXPECT_EQ((uint32_t)((N + 63) / 64), host_grid[0]);
    EXPECT_EQ(1u, host_grid[1]);
    EXPECT_EQ(1u, host_grid[2]);
    for (int i = 0; i < N; ++i) EXPECT_EQ((uint32_t)(i + 1000), host[i]);

    md_gpu_free(data,  f.compute);
    md_gpu_free(ones,  f.compute);
    md_gpu_free(count, f.compute);
    md_gpu_free(grid,  f.compute);
    gpu_close(&f);
}

/* =========================================================================
   Graphs
   ========================================================================= */

UTEST(gpu, graph_capture_and_replay) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 1024 };
    uint32_t* a = (uint32_t*)md_gpu_malloc(f.pool, N * sizeof(uint32_t), f.compute);
    uint32_t* b = (uint32_t*)md_gpu_malloc(f.pool, N * sizeof(uint32_t), f.compute);
    ASSERT_TRUE(a && b);

    ASSERT_TRUE(md_gpu_capture_begin(f.compute, "test graph"));
    ASSERT_TRUE(md_gpu_is_capturing(f.compute));

    EXPECT_EQ(0u, md_gpu_capture_next_index(f.compute));
    fill_args_t fa = {0};
    fa.n = N; fa.base = 0; fa.dst = DEV_ADDR(a);
    ASSERT_TRUE(md_gpu_launch(f.compute, f.k_fill, md_gpu_grid_1d(N, 64), &fa, sizeof(fa)));

    uint32_t scale_index = md_gpu_capture_next_index(f.compute);
    scale_args_t sa = {0};
    sa.n = N; sa.mul = 1; sa.add = 7;
    sa.src = DEV_ADDR(a); sa.dst = DEV_ADDR(b);
    ASSERT_TRUE(md_gpu_launch(f.compute, f.k_scale, md_gpu_grid_1d(N, 64), &sa, sizeof(sa)));

    md_gpu_graph_t g = md_gpu_capture_end(f.compute);
    ASSERT_TRUE(g != NULL);
    ASSERT_FALSE(md_gpu_is_capturing(f.compute));
    EXPECT_EQ(2u, md_gpu_graph_launch_count(g));

    /* First launch. */
    ASSERT_TRUE(md_gpu_graph_launch(g, f.compute));
    uint32_t host[N];
    ASSERT_TRUE(md_gpu_memcpy_async(host, b, sizeof(host), f.compute));
    md_gpu_stream_sync(f.compute);
    for (int i = 0; i < N; ++i) EXPECT_EQ((uint32_t)(i + 7), host[i]);

    /* Patch the argument block in place and replay -- no re-recording. */
    scale_args_t* patched = (scale_args_t*)md_gpu_graph_args(g, scale_index);
    ASSERT_TRUE(patched != NULL);
    EXPECT_EQ(7u, patched->add);
    patched->add = 1000;
    patched->mul = 2;

    ASSERT_TRUE(md_gpu_graph_launch(g, f.compute));
    memset(host, 0, sizeof(host));
    ASSERT_TRUE(md_gpu_memcpy_async(host, b, sizeof(host), f.compute));
    md_gpu_stream_sync(f.compute);
    for (int i = 0; i < N; ++i) EXPECT_EQ((uint32_t)(i * 2 + 1000), host[i]);

    md_gpu_graph_destroy(g);
    md_gpu_free(a, f.compute);
    md_gpu_free(b, f.compute);
    gpu_close(&f);
}

UTEST(gpu, graph_replayed_many_times) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 256, REPLAYS = 16 };
    uint32_t* d = (uint32_t*)md_gpu_malloc(f.pool, N * sizeof(uint32_t), f.compute);
    ASSERT_TRUE(d != NULL);

    ASSERT_TRUE(md_gpu_memset_async(d, 0, N * sizeof(uint32_t), f.compute));
    md_gpu_stream_sync(f.compute);

    ASSERT_TRUE(md_gpu_capture_begin(f.compute, "increment"));
    bump_args_t ba = {0};
    ba.n = N; ba.delta = 1; ba.dst = DEV_ADDR(d);
    ASSERT_TRUE(md_gpu_launch(f.compute, f.k_bump, md_gpu_grid_1d(N, 64), &ba, sizeof(ba)));
    md_gpu_graph_t g = md_gpu_capture_end(f.compute);
    ASSERT_TRUE(g != NULL);

    for (int i = 0; i < REPLAYS; ++i) {
        ASSERT_TRUE(md_gpu_graph_launch(g, f.compute));
    }

    uint32_t host[N];
    ASSERT_TRUE(md_gpu_memcpy_async(host, d, sizeof(host), f.compute));
    md_gpu_stream_sync(f.compute);
    for (int i = 0; i < N; ++i) EXPECT_EQ((uint32_t)REPLAYS, host[i]);

    md_gpu_graph_destroy(g);
    md_gpu_free(d, f.compute);
    gpu_close(&f);
}

/* =========================================================================
   Host callbacks
   ========================================================================= */

typedef struct {
    int      fired;
    uint32_t observed;
    uint32_t* watch;
} cb_state_t;

static void gpu_test_callback(void* user) {
    cb_state_t* st = (cb_state_t*)user;
    st->fired++;
    st->observed = st->watch ? st->watch[0] : 0;
}

UTEST(gpu, host_callback_fires_in_poll) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 512 };
    uint32_t* d = (uint32_t*)md_gpu_malloc(f.pool, N * sizeof(uint32_t), f.compute);
    ASSERT_TRUE(d != NULL);

    static uint32_t readback[N];
    memset(readback, 0, sizeof(readback));

    cb_state_t st = {0};
    st.watch = readback;

    fill_args_t a = {0};
    a.n = N; a.base = 77; a.dst = DEV_ADDR(d);
    ASSERT_TRUE(md_gpu_launch(f.compute, f.k_fill, md_gpu_grid_1d(N, 64), &a, sizeof(a)));
    ASSERT_TRUE(md_gpu_memcpy_async(readback, d, sizeof(readback), f.compute));

    /* Nothing has run the callback yet -- it only fires inside device_poll. */
    ASSERT_TRUE(md_gpu_launch_host_fn(f.compute, gpu_test_callback, &st));
    EXPECT_EQ(0, st.fired);

    md_gpu_stream_sync(f.compute);
    EXPECT_EQ(0, st.fired);       /* still not fired: sync is not poll */

    uint32_t fired = md_gpu_device_poll(f.dev);
    EXPECT_GT(fired, 0u);
    EXPECT_EQ(1, st.fired);
    /* The device-to-host copy issued before the callback is visible to it. */
    EXPECT_EQ(77u, st.observed);

    /* Polling again does not re-fire. */
    md_gpu_device_poll(f.dev);
    EXPECT_EQ(1, st.fired);

    md_gpu_free(d, f.compute);
    gpu_close(&f);
}

/* A host callback must observe every write issued into the stream before it,
   including a device-to-host copy that md_gpu completes with an internal
   callback of its own.

   This is a race, not a logic error: the copy's completion callback and the
   user callback are registered against the same sync value, and readiness used
   to be decided by re-reading the timeline per entry. If the timeline advanced
   between the two reads, the user callback overtook the copy and saw stale
   host memory. It reproduced roughly one run in four, so the loop matters. */
typedef struct {
    const uint32_t* watch;
    uint32_t        observed;
    int             fired;
} order_state_t;

static void gpu_order_callback(void* user) {
    order_state_t* st = (order_state_t*)user;
    st->observed = st->watch[0];
    st->fired    = 1;
}

UTEST(gpu, host_callback_observes_preceding_readback) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 256, ITERATIONS = 32 };
    uint32_t* d = (uint32_t*)md_gpu_malloc(f.pool, N * sizeof(uint32_t), f.compute);
    ASSERT_TRUE(d != NULL);

    static uint32_t host[N];

    for (uint32_t iter = 0; iter < ITERATIONS; ++iter) {
        memset(host, 0, sizeof(host));

        fill_args_t a = {0};
        a.n = N; a.base = 1000 + iter; a.dst = DEV_ADDR(d);
        ASSERT_TRUE(md_gpu_launch(f.compute, f.k_fill, md_gpu_grid_1d(N, 64), &a, sizeof(a)));

        order_state_t st = {0};
        st.watch = host;

        ASSERT_TRUE(md_gpu_memcpy_async(host, d, sizeof(host), f.compute));
        ASSERT_TRUE(md_gpu_launch_host_fn(f.compute, gpu_order_callback, &st));

        while (!st.fired) md_gpu_device_poll(f.dev);

        /* The callback ran after the copy, so it saw the fill's output. */
        EXPECT_EQ(1000u + iter, st.observed);
    }

    md_gpu_free(d, f.compute);
    gpu_close(&f);
}

UTEST(gpu, sync_on_complete) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 64 };
    uint32_t* d = (uint32_t*)md_gpu_malloc(f.pool, N * sizeof(uint32_t), f.compute);
    ASSERT_TRUE(d != NULL);

    fill_args_t a = {0};
    a.n = N; a.base = 1; a.dst = DEV_ADDR(d);
    ASSERT_TRUE(md_gpu_launch(f.compute, f.k_fill, md_gpu_grid_1d(N, 64), &a, sizeof(a)));

    md_gpu_sync_t sync = md_gpu_stream_record(f.compute);
    cb_state_t st = {0};
    ASSERT_TRUE(md_gpu_sync_on_complete(sync, gpu_test_callback, &st));

    md_gpu_sync_wait(sync);
    EXPECT_TRUE(md_gpu_sync_is_complete(sync));
    md_gpu_device_poll(f.dev);
    EXPECT_EQ(1, st.fired);

    md_gpu_free(d, f.compute);
    gpu_close(&f);
}

/* =========================================================================
   Kernel introspection
   ========================================================================= */

UTEST(gpu, kernel_info_recovers_group_size) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_kernel_info_t info = {0};
    ASSERT_TRUE(md_gpu_kernel_info(f.k_fill, &info));
    EXPECT_EQ(64u, info.group_size[0]);
    EXPECT_EQ(1u,  info.group_size[1]);
    EXPECT_EQ(1u,  info.group_size[2]);
    EXPECT_GT(info.max_threads_per_group, 0u);

    ASSERT_TRUE(md_gpu_kernel_info(f.k_tex_write, &info));
    EXPECT_EQ(4u, info.group_size[0]);
    EXPECT_EQ(4u, info.group_size[1]);
    EXPECT_EQ(4u, info.group_size[2]);

    gpu_close(&f);
}

/* =========================================================================
   Stress: many allocations and launches, checking for leaks and wall-hits
   ========================================================================= */

UTEST(gpu, stress_alloc_launch_free_cycles) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 4096, CYCLES = 64 };
    md_gpu_pool_t pool = md_gpu_pool_create(f.dev, &(md_gpu_pool_desc_t){ .flags = MD_GPU_MEM_DEVICE, .release_threshold = 32 * 1024 * 1024, .label = "stress" });
    ASSERT_TRUE(pool != NULL);

    for (int c = 0; c < CYCLES; ++c) {
        uint32_t* d = (uint32_t*)md_gpu_malloc(pool, N * sizeof(uint32_t), f.compute);
        ASSERT_TRUE(d != NULL);
        fill_args_t a = {0};
        a.n = N; a.base = (uint32_t)c; a.dst = DEV_ADDR(d);
        ASSERT_TRUE(md_gpu_launch(f.compute, f.k_fill, md_gpu_grid_1d(N, 64), &a, sizeof(a)));
        md_gpu_free(d, f.compute);
        if ((c & 7) == 0) md_gpu_device_poll(f.dev);
    }
    md_gpu_stream_sync(f.compute);

    /* Every block was freed on the same stream, so the pool should have
       recycled one block rather than committing 64. */
    md_gpu_pool_stats_t st = {0};
    md_gpu_pool_stats(pool, &st);
    EXPECT_EQ(0u, (unsigned)st.bytes_in_use);
    EXPECT_LE(st.bytes_reserved, (uint64_t)(N * sizeof(uint32_t) * 4));
    /* Almost every allocation should have been served from cache. */
    EXPECT_GE(st.reuse_count, st.alloc_count - 1);

    md_gpu_pool_destroy(pool);
    gpu_close(&f);
}

/* Concurrency comes from streams, not from a flag. Two streams writing
   disjoint buffers must both land, and each stream is internally ordered. */
UTEST(gpu, two_streams_are_independent) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 512 };
    md_gpu_stream_t sa = md_gpu_stream_create(f.dev, MD_GPU_STREAM_COMPUTE, "a");
    md_gpu_stream_t sb = md_gpu_stream_create(f.dev, MD_GPU_STREAM_COMPUTE, "b");
    ASSERT_TRUE(sa && sb);

    uint32_t* a = (uint32_t*)md_gpu_malloc(f.pool, N * sizeof(uint32_t), sa);
    uint32_t* b = (uint32_t*)md_gpu_malloc(f.pool, N * sizeof(uint32_t), sb);
    ASSERT_TRUE(a && b);

    fill_args_t fa = {0}; fa.n = N; fa.base = 10; fa.dst = DEV_ADDR(a);
    fill_args_t fb = {0}; fb.n = N; fb.base = 20; fb.dst = DEV_ADDR(b);
    ASSERT_TRUE(md_gpu_launch(sa, f.k_fill, md_gpu_grid_1d(N, 64), &fa, sizeof(fa)));
    ASSERT_TRUE(md_gpu_launch(sb, f.k_fill, md_gpu_grid_1d(N, 64), &fb, sizeof(fb)));

    /* Each stream doubles its own buffer in place -- ordered within a stream. */
    scale_args_t da = {0}; da.n = N; da.mul = 2; da.src = DEV_ADDR(a); da.dst = DEV_ADDR(a);
    scale_args_t db = {0}; db.n = N; db.mul = 3; db.src = DEV_ADDR(b); db.dst = DEV_ADDR(b);
    ASSERT_TRUE(md_gpu_launch(sa, f.k_scale, md_gpu_grid_1d(N, 64), &da, sizeof(da)));
    ASSERT_TRUE(md_gpu_launch(sb, f.k_scale, md_gpu_grid_1d(N, 64), &db, sizeof(db)));

    uint32_t ha[N], hb[N];
    ASSERT_TRUE(md_gpu_memcpy_async(ha, a, sizeof(ha), sa));
    ASSERT_TRUE(md_gpu_memcpy_async(hb, b, sizeof(hb), sb));
    md_gpu_stream_sync(sa);
    md_gpu_stream_sync(sb);

    for (int i = 0; i < N; ++i) {
        EXPECT_EQ((uint32_t)((10 + i) * 2), ha[i]);
        EXPECT_EQ((uint32_t)((20 + i) * 3), hb[i]);
    }

    md_gpu_free(a, sa);
    md_gpu_free(b, sb);
    md_gpu_stream_destroy(sa);
    md_gpu_stream_destroy(sb);
    gpu_close(&f);
}

/* =========================================================================
   Argument-struct layout

   md_gpu.h spends a whole section on md_gpu_float4 and friends because SPIR-V
   and MSL lay vectors out differently, and nothing exercised them. A mismatch
   here does not announce itself -- it shows up as a wrong number several
   launches downstream -- so read the fields straight back.
   ========================================================================= */

UTEST(gpu, arg_struct_layout_matches_shader) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 10 };
    uint32_t* d = (uint32_t*)md_gpu_malloc(f.pool, N * sizeof(uint32_t), f.compute);
    ASSERT_TRUE(d != NULL);
    ASSERT_TRUE(md_gpu_memset_async(d, 0xEE, N * sizeof(uint32_t), f.compute));

    layout_probe_args_t a = {0};
    a.dim_x = 11; a.dim_y = 22; a.dim_z = 33;
    a.scale = 1.5f;
    a.pair.x = 44; a.pair.y = 55;
    a.v4.x = 2.5f; a.v4.y = 3.5f; a.v4.z = 4.5f; a.v4.w = 5.5f;
    a.dst = DEV_ADDR(d);
    ASSERT_TRUE(md_gpu_launch(f.compute, f.k_layout, md_gpu_grid(1, 1, 1), &a, sizeof(a)));

    uint32_t host[N];
    ASSERT_TRUE(md_gpu_memcpy_async(host, d, sizeof(host), f.compute));
    md_gpu_stream_sync(f.compute);

    float scale, v4[4];
    memcpy(&scale, &host[3], sizeof(scale));
    memcpy(v4, &host[6], sizeof(v4));

    EXPECT_EQ(11u, host[0]);
    EXPECT_EQ(22u, host[1]);
    EXPECT_EQ(33u, host[2]);
    EXPECT_EQ(1.5f, scale);          /* scalar tail after the 8-aligned uint2 */
    EXPECT_EQ(44u, host[4]);
    EXPECT_EQ(55u, host[5]);         /* uint2: same size, different alignment */
    EXPECT_EQ(2.5f, v4[0]);
    EXPECT_EQ(3.5f, v4[1]);
    EXPECT_EQ(4.5f, v4[2]);
    EXPECT_EQ(5.5f, v4[3]);

    md_gpu_free(d, f.compute);
    gpu_close(&f);
}

/* A texture handle lives inside the argument struct, and the two backends
   resolve it by entirely different means. This separates the two things that
   can go wrong: whether the struct survived the handle (dst[0]) and whether the
   handle itself resolved to a real texture (the readback). */
UTEST(gpu, texture_handle_reaches_shader) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { D = 4, VOXELS = D * D * D, MARKER = 0xC0FFEEu };
    md_gpu_tex_desc_t td = {0};
    td.width = D; td.height = D; td.depth = D;
    td.format = MD_GPU_FORMAT_R32_FLOAT;
    td.flags  = MD_GPU_TEX_STORAGE;
    td.label  = "probe";
    md_gpu_tex_t tex = md_gpu_tex_create(f.dev, &td);
    ASSERT_TRUE(tex != 0);

    uint32_t* d = (uint32_t*)md_gpu_malloc(f.pool, 2 * sizeof(uint32_t), f.compute);
    ASSERT_TRUE(d != NULL);
    ASSERT_TRUE(md_gpu_memset_async(d, 0, 2 * sizeof(uint32_t), f.compute));

    float zero[VOXELS] = {0};
    ASSERT_TRUE(md_gpu_memcpy_to_tex_async(tex, NULL, zero, sizeof(zero), f.compute));

    tex_probe_args_t a = {0};
    a.n      = 7;
    a.dst    = DEV_ADDR(d);
    a.tex    = tex;
    a.marker = MARKER;
    ASSERT_TRUE(md_gpu_launch(f.compute, f.k_tex_probe, md_gpu_grid(1, 1, 1), &a, sizeof(a)));

    uint32_t host[2] = {0};
    float voxels[VOXELS];
    memset(voxels, 0xFF, sizeof(voxels));
    ASSERT_TRUE(md_gpu_memcpy_async(host, d, sizeof(host), f.compute));
    ASSERT_TRUE(md_gpu_memcpy_from_tex_async(voxels, tex, NULL, sizeof(voxels), f.compute));
    md_gpu_stream_sync(f.compute);

    /* Members after the handle are still where the C mirror put them. */
    EXPECT_EQ((uint32_t)MARKER, host[0]);
    EXPECT_EQ(7u, host[1]);
    /* And the handle resolved: the shader wrote through it. */
    EXPECT_EQ(1.0f, voxels[0]);
    EXPECT_EQ(2.0f, voxels[1]);

    md_gpu_free(d, f.compute);
    md_gpu_tex_destroy(tex, f.compute);
    md_gpu_device_poll(f.dev);
    gpu_close(&f);
}

/* =========================================================================
   Staging and the transient arena
   ========================================================================= */

/* Bigger than the backends' transient page size, so upload_begin has to commit
   a fresh page rather than carve one up. */
UTEST(gpu, upload_larger_than_one_arena_page) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 512 * 1024 };          /* 2 MiB, well past a 256 KiB page */
    uint32_t* d = (uint32_t*)md_gpu_malloc(f.pool, N * sizeof(uint32_t), f.compute);
    ASSERT_TRUE(d != NULL);

    uint32_t* p = (uint32_t*)md_gpu_upload_begin(f.compute, d, N * sizeof(uint32_t));
    ASSERT_TRUE(p != NULL);
    for (int i = 0; i < N; ++i) p[i] = (uint32_t)i * 2654435761u;
    ASSERT_TRUE(md_gpu_upload_end(f.compute));

    static uint32_t host[N];
    ASSERT_TRUE(md_gpu_memcpy_async(host, d, sizeof(host), f.compute));
    md_gpu_stream_sync(f.compute);
    for (int i = 0; i < N; ++i) ASSERT_EQ((uint32_t)i * 2654435761u, host[i]);

    md_gpu_free(d, f.compute);
    gpu_close(&f);
}

/* A device-to-host copy of non-pinned memory lands in transient staging and is
   delivered later by an internal callback. Issuing more work -- which draws
   from that same arena -- before the delivery must not recycle the page out
   from under it. The interleaving is what makes this bite, so alternate. */
UTEST(gpu, readback_survives_subsequent_staged_work) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 256, ROUNDS = 24 };
    uint32_t* a = (uint32_t*)md_gpu_malloc(f.pool, N * sizeof(uint32_t), f.compute);
    uint32_t* b = (uint32_t*)md_gpu_malloc(f.pool, N * sizeof(uint32_t), f.compute);
    ASSERT_TRUE(a && b);

    static uint32_t upload[N];
    static uint32_t host[N];

    for (uint32_t r = 0; r < ROUNDS; ++r) {
        fill_args_t fa = {0};
        fa.n = N; fa.base = 3000 + r; fa.dst = DEV_ADDR(a);
        ASSERT_TRUE(md_gpu_launch(f.compute, f.k_fill, md_gpu_grid_1d(N, 64), &fa, sizeof(fa)));

        memset(host, 0, sizeof(host));
        ASSERT_TRUE(md_gpu_memcpy_async(host, a, sizeof(host), f.compute));

        /* Host-to-device staging issued *after* the readback, before it has
           been delivered. Same arena, same stream. */
        for (int i = 0; i < N; ++i) upload[i] = 0xDEADBEEFu;
        ASSERT_TRUE(md_gpu_memcpy_async(b, upload, sizeof(upload), f.compute));

        md_gpu_stream_sync(f.compute);
        md_gpu_device_poll(f.dev);

        for (int i = 0; i < N; ++i) ASSERT_EQ((uint32_t)(3000 + r + i), host[i]);
    }

    md_gpu_free(a, f.compute);
    md_gpu_free(b, f.compute);
    gpu_close(&f);
}

/* =========================================================================
   Graphs beyond dispatches
   ========================================================================= */

/* A graph is replayed many times and long after capture. Anything it recorded
   that points into per-stream transient storage will have been recycled by
   then, so a captured upload has to own its bytes. */
UTEST(gpu, graph_replays_captured_copies_and_fills) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 256, REPLAYS = 8 };
    uint32_t* d = (uint32_t*)md_gpu_malloc(f.pool, N * sizeof(uint32_t), f.compute);
    uint32_t* e = (uint32_t*)md_gpu_malloc(f.pool, N * sizeof(uint32_t), f.compute);
    ASSERT_TRUE(d && e);

    uint32_t src[N];
    for (int i = 0; i < N; ++i) src[i] = (uint32_t)(i + 1);

    ASSERT_TRUE(md_gpu_capture_begin(f.compute, "copies"));
    ASSERT_TRUE(md_gpu_memset_async(d, 0, N * sizeof(uint32_t), f.compute));
    ASSERT_TRUE(md_gpu_memcpy_async(d, src, sizeof(src), f.compute));   /* host -> device */
    ASSERT_TRUE(md_gpu_memcpy_async(e, d, N * sizeof(uint32_t), f.compute));  /* device -> device */
    md_gpu_graph_t g = md_gpu_capture_end(f.compute);
    ASSERT_TRUE(g != NULL);

    /* Churn the stream's transient arena before the first replay: if the graph
       recorded a pointer into it, this is what corrupts it. */
    memset(src, 0, sizeof(src));
    for (int i = 0; i < 8; ++i) {
        ASSERT_TRUE(md_gpu_memcpy_async(e, src, sizeof(src), f.compute));
    }
    md_gpu_stream_sync(f.compute);
    for (int i = 0; i < N; ++i) src[i] = (uint32_t)(i + 1);   /* restore host copy only */

    uint32_t host[N];
    for (int r = 0; r < REPLAYS; ++r) {
        ASSERT_TRUE(md_gpu_graph_launch(g, f.compute));
        memset(host, 0, sizeof(host));
        ASSERT_TRUE(md_gpu_memcpy_async(host, e, sizeof(host), f.compute));
        md_gpu_stream_sync(f.compute);
        for (int i = 0; i < N; ++i) ASSERT_EQ((uint32_t)(i + 1), host[i]);
    }

    md_gpu_graph_destroy(g);
    md_gpu_free(d, f.compute);
    md_gpu_free(e, f.compute);
    gpu_close(&f);
}

/* =========================================================================
   Texture regions and samplers
   ========================================================================= */

/* Every texture copy so far has covered the whole image, so the region path --
   origin, extent, row and slice pitch -- is entirely unexercised. */
UTEST(gpu, texture_region_subrange_roundtrip) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { D = 8, SUB = 4, VOXELS = D * D * D, SUBVOX = SUB * SUB * SUB };
    md_gpu_tex_desc_t td = {0};
    td.width = D; td.height = D; td.depth = D;
    td.format = MD_GPU_FORMAT_R32_FLOAT;
    td.flags  = MD_GPU_TEX_STORAGE;
    md_gpu_tex_t tex = md_gpu_tex_create(f.dev, &td);
    ASSERT_TRUE(tex != 0);

    float base[VOXELS];
    for (int i = 0; i < VOXELS; ++i) base[i] = (float)i;
    ASSERT_TRUE(md_gpu_memcpy_to_tex_async(tex, NULL, base, sizeof(base), f.compute));

    /* Overwrite one corner octant. */
    float patch[SUBVOX];
    for (int i = 0; i < SUBVOX; ++i) patch[i] = -1.0f - (float)i;
    md_gpu_tex_region_t region = {0};
    region.offset[0] = SUB; region.offset[1] = SUB; region.offset[2] = SUB;
    region.extent[0] = SUB; region.extent[1] = SUB; region.extent[2] = SUB;
    ASSERT_TRUE(md_gpu_memcpy_to_tex_async(tex, &region, patch, sizeof(patch), f.compute));

    /* Read the same octant back and confirm it is what we wrote. */
    float back[SUBVOX];
    memset(back, 0, sizeof(back));
    ASSERT_TRUE(md_gpu_memcpy_from_tex_async(back, tex, &region, sizeof(back), f.compute));

    /* And the rest of the volume is untouched: sample the opposite corner. */
    md_gpu_tex_region_t other = {0};
    other.extent[0] = SUB; other.extent[1] = SUB; other.extent[2] = SUB;
    float corner[SUBVOX];
    memset(corner, 0, sizeof(corner));
    ASSERT_TRUE(md_gpu_memcpy_from_tex_async(corner, tex, &other, sizeof(corner), f.compute));
    md_gpu_stream_sync(f.compute);

    for (int i = 0; i < SUBVOX; ++i) EXPECT_EQ(-1.0f - (float)i, back[i]);
    for (int z = 0; z < SUB; ++z) {
        for (int y = 0; y < SUB; ++y) {
            for (int x = 0; x < SUB; ++x) {
                float expect = (float)(z * D * D + y * D + x);
                EXPECT_EQ(expect, corner[z * SUB * SUB + y * SUB + x]);
            }
        }
    }

    md_gpu_tex_destroy(tex, f.compute);
    md_gpu_device_poll(f.dev);
    gpu_close(&f);
}

/* =========================================================================
   Execution hazards

   md_gpu.h's central promise is that every operation observes all writes made
   by the operations before it in the same stream. Neither backend gets that for
   free, and both lose it in a way no small test notices.

   Metal decides whether two encoders inside a command buffer may overlap from
   hazard tracking, and this API declares nothing per dispatch -- buffers arrive
   as raw device addresses, textures as resource ids inside the argument struct,
   and MTLResidencySet supplies residency and explicitly not usage. So the
   dependency graph is empty, adjacent encoders look independent, and a readback
   blit is free to start copying before the dispatch that fills the volume has
   finished. Vulkan has the same shape of problem answered differently: the
   ordering is whatever vkCmdPipelineBarrier the backend emits, and a missing or
   over-narrow stage/access mask fails exactly here.

   The two rules these tests exist to enforce:

     - SIZE. A producer must run long enough for an unordered consumer to get
       ahead of it. The 8^3 volume in the texture tests above completes before a
       blit could possibly overtake it, which is why this bug survived a passing
       suite and turned up as artifacts in a real GTO volume instead.

     - ROUNDS. Each round writes a different value, so a stale read reports the
       *previous round's* value rather than zero or garbage. "Expected 3000, got
       2000" says ordering; "expected 3000, got 0" says something else entirely.

   The matrix below covers producer x consumer x separation. Every case is
   backend-agnostic: nothing here mentions fences, barriers or encoders.
   ========================================================================= */

enum {
    HAZ_D      = 128,                        /* 128^3 volume, 8 MiB          */
    HAZ_VOX    = HAZ_D * HAZ_D * HAZ_D,
    HAZ_N      = 1 << 20,                    /* 1 M element buffer, 4 MiB    */
    HAZ_ROUNDS = 4,
    HAZ_STRIDE = 389,                        /* coprime with everything here */
};

/* dispatch -> blit, buffer. The dispatch writes through a device address the
   backend never declared; the blit reads the same range. */
UTEST(gpu, hazard_dispatch_to_blit_buffer) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    uint32_t* host = (uint32_t*)malloc(HAZ_N * sizeof(uint32_t));
    ASSERT_TRUE(host != NULL);
    uint32_t* d = (uint32_t*)md_gpu_malloc(f.pool, HAZ_N * sizeof(uint32_t), f.compute);
    ASSERT_TRUE(d != NULL);

    for (int r = 0; r < HAZ_ROUNDS; ++r) {
        const uint32_t base = (uint32_t)(r + 1) * 1000000u;
        fill_args_t fa = {0};
        fa.n = HAZ_N; fa.base = base; fa.dst = DEV_ADDR(d);
        ASSERT_TRUE(md_gpu_launch(f.compute, f.k_fill, md_gpu_grid_1d(HAZ_N, 64), &fa, sizeof(fa)));

        /* No flush: producer and consumer land in one command buffer, which is
           where the encoders may overlap. */
        memset(host, 0, HAZ_N * sizeof(uint32_t));
        ASSERT_TRUE(md_gpu_memcpy_async(host, d, HAZ_N * sizeof(uint32_t), f.compute));
        md_gpu_stream_sync(f.compute);
        md_gpu_device_poll(f.dev);

        for (int i = 0; i < HAZ_N; i += HAZ_STRIDE) ASSERT_EQ(base + (uint32_t)i, host[i]);
    }

    md_gpu_free(d, f.compute);
    free(host);
    gpu_close(&f);
}

/* dispatch -> blit, texture. The GTO volume readback, reduced. */
UTEST(gpu, hazard_dispatch_to_blit_texture) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_tex_desc_t td = {0};
    td.width = HAZ_D; td.height = HAZ_D; td.depth = HAZ_D;
    td.format = MD_GPU_FORMAT_R32_FLOAT;
    td.flags  = MD_GPU_TEX_STORAGE;
    td.label  = "hazard volume";
    md_gpu_tex_t tex = md_gpu_tex_create(f.dev, &td);
    ASSERT_TRUE(tex != 0);

    float* host = (float*)malloc(HAZ_VOX * sizeof(float));
    ASSERT_TRUE(host != NULL);

    for (int r = 0; r < HAZ_ROUNDS; ++r) {
        const float scale = (float)(r + 1);
        tex_args_t ta = {0};
        ta.dim[0] = HAZ_D; ta.dim[1] = HAZ_D; ta.dim[2] = HAZ_D;
        ta.tex = tex; ta.scale = scale;
        ASSERT_TRUE(md_gpu_launch(f.compute, f.k_tex_write,
                                  md_gpu_grid(HAZ_D/4, HAZ_D/4, HAZ_D/4), &ta, sizeof(ta)));

        memset(host, 0, HAZ_VOX * sizeof(float));
        ASSERT_TRUE(md_gpu_memcpy_from_tex_async(host, tex, NULL, HAZ_VOX * sizeof(float), f.compute));
        md_gpu_stream_sync(f.compute);
        md_gpu_device_poll(f.dev);

        for (int i = 0; i < HAZ_VOX; i += HAZ_STRIDE) ASSERT_EQ((float)i * scale, host[i]);
    }

    free(host);
    md_gpu_tex_destroy(tex, f.compute);
    md_gpu_device_poll(f.dev);
    gpu_close(&f);
}

/* blit -> dispatch, buffer. The upload must be complete before the kernel that
   consumes it starts. */
UTEST(gpu, hazard_blit_to_dispatch_buffer) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    uint32_t* src  = (uint32_t*)malloc(HAZ_N * sizeof(uint32_t));
    uint32_t* host = (uint32_t*)malloc(HAZ_N * sizeof(uint32_t));
    ASSERT_TRUE(src && host);
    uint32_t* a = (uint32_t*)md_gpu_malloc(f.pool, HAZ_N * sizeof(uint32_t), f.compute);
    uint32_t* b = (uint32_t*)md_gpu_malloc(f.pool, HAZ_N * sizeof(uint32_t), f.compute);
    ASSERT_TRUE(a && b);

    for (int r = 0; r < HAZ_ROUNDS; ++r) {
        const uint32_t tag = (uint32_t)(r + 1) * 7u;
        for (int i = 0; i < HAZ_N; ++i) src[i] = (uint32_t)i + tag;

        ASSERT_TRUE(md_gpu_memcpy_async(a, src, HAZ_N * sizeof(uint32_t), f.compute));

        scale_args_t sa = {0};
        sa.n = HAZ_N; sa.mul = 1; sa.add = 0;
        sa.src = DEV_ADDR(a); sa.dst = DEV_ADDR(b);
        ASSERT_TRUE(md_gpu_launch(f.compute, f.k_scale, md_gpu_grid_1d(HAZ_N, 64), &sa, sizeof(sa)));

        memset(host, 0, HAZ_N * sizeof(uint32_t));
        ASSERT_TRUE(md_gpu_memcpy_async(host, b, HAZ_N * sizeof(uint32_t), f.compute));
        md_gpu_stream_sync(f.compute);
        md_gpu_device_poll(f.dev);

        for (int i = 0; i < HAZ_N; i += HAZ_STRIDE) ASSERT_EQ((uint32_t)i + tag, host[i]);
    }

    md_gpu_free(a, f.compute);
    md_gpu_free(b, f.compute);
    free(src); free(host);
    gpu_close(&f);
}

/* blit -> dispatch, texture. */
UTEST(gpu, hazard_blit_to_dispatch_texture) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { D = 64, VOX = D * D * D };
    md_gpu_tex_desc_t td = {0};
    td.width = D; td.height = D; td.depth = D;
    td.format = MD_GPU_FORMAT_R32_FLOAT;
    td.flags  = MD_GPU_TEX_STORAGE;
    md_gpu_tex_t tex = md_gpu_tex_create(f.dev, &td);
    ASSERT_TRUE(tex != 0);

    float* src  = (float*)malloc(VOX * sizeof(float));
    float* host = (float*)malloc(VOX * sizeof(float));
    ASSERT_TRUE(src && host);
    float* d = (float*)md_gpu_malloc(f.pool, VOX * sizeof(float), f.compute);
    ASSERT_TRUE(d != NULL);

    for (int r = 0; r < HAZ_ROUNDS; ++r) {
        const float tag = (float)(r + 1) * 1000.0f;
        for (int i = 0; i < VOX; ++i) src[i] = (float)i + tag;

        ASSERT_TRUE(md_gpu_memcpy_to_tex_async(tex, NULL, src, VOX * sizeof(float), f.compute));

        tex_read_args_t ra = {0};
        ra.dim[0] = D; ra.dim[1] = D; ra.dim[2] = D;
        ra.tex = tex; ra.dst = DEV_ADDR(d);
        ASSERT_TRUE(md_gpu_launch(f.compute, f.k_tex_read, md_gpu_grid(D/4, D/4, D/4), &ra, sizeof(ra)));

        memset(host, 0, VOX * sizeof(float));
        ASSERT_TRUE(md_gpu_memcpy_async(host, d, VOX * sizeof(float), f.compute));
        md_gpu_stream_sync(f.compute);
        md_gpu_device_poll(f.dev);

        for (int i = 0; i < VOX; i += 61) ASSERT_EQ(src[i], host[i]);
    }

    md_gpu_free(d, f.compute);
    free(src); free(host);
    md_gpu_tex_destroy(tex, f.compute);
    md_gpu_device_poll(f.dev);
    gpu_close(&f);
}

/* memset -> dispatch. A fill is a transfer operation like any other, and a
   kernel that accumulates into the cleared range must not start early. */
UTEST(gpu, hazard_memset_to_dispatch) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    uint32_t* host = (uint32_t*)malloc(HAZ_N * sizeof(uint32_t));
    ASSERT_TRUE(host != NULL);
    uint32_t* d = (uint32_t*)md_gpu_malloc(f.pool, HAZ_N * sizeof(uint32_t), f.compute);
    ASSERT_TRUE(d != NULL);

    for (int r = 0; r < HAZ_ROUNDS; ++r) {
        const uint8_t  byte  = (uint8_t)(r + 1);
        const uint32_t word  = (uint32_t)byte * 0x01010101u;
        const uint32_t delta = (uint32_t)(r + 1) * 13u;

        ASSERT_TRUE(md_gpu_memset_async(d, byte, HAZ_N * sizeof(uint32_t), f.compute));

        bump_args_t ba = {0};
        ba.n = HAZ_N; ba.delta = delta; ba.dst = DEV_ADDR(d);
        ASSERT_TRUE(md_gpu_launch(f.compute, f.k_bump, md_gpu_grid_1d(HAZ_N, 64), &ba, sizeof(ba)));

        memset(host, 0, HAZ_N * sizeof(uint32_t));
        ASSERT_TRUE(md_gpu_memcpy_async(host, d, HAZ_N * sizeof(uint32_t), f.compute));
        md_gpu_stream_sync(f.compute);
        md_gpu_device_poll(f.dev);

        for (int i = 0; i < HAZ_N; i += HAZ_STRIDE) ASSERT_EQ(word + delta, host[i]);
    }

    md_gpu_free(d, f.compute);
    free(host);
    gpu_close(&f);
}

/* dispatch -> dispatch through a texture. Both sides are compute, so on Metal
   this is the serial-encoder path rather than the encoder boundary -- a
   different mechanism reaching the same guarantee, and worth pinning
   separately. */
UTEST(gpu, hazard_dispatch_to_dispatch_via_texture) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { D = 64, VOX = D * D * D };
    md_gpu_tex_desc_t td = {0};
    td.width = D; td.height = D; td.depth = D;
    td.format = MD_GPU_FORMAT_R32_FLOAT;
    td.flags  = MD_GPU_TEX_STORAGE;
    md_gpu_tex_t tex = md_gpu_tex_create(f.dev, &td);
    ASSERT_TRUE(tex != 0);

    float* host = (float*)malloc(VOX * sizeof(float));
    ASSERT_TRUE(host != NULL);
    float* d = (float*)md_gpu_malloc(f.pool, VOX * sizeof(float), f.compute);
    ASSERT_TRUE(d != NULL);

    for (int r = 0; r < HAZ_ROUNDS; ++r) {
        const float scale = (float)(r + 1) * 0.5f;

        tex_args_t ta = {0};
        ta.dim[0] = D; ta.dim[1] = D; ta.dim[2] = D;
        ta.tex = tex; ta.scale = scale;
        ASSERT_TRUE(md_gpu_launch(f.compute, f.k_tex_write, md_gpu_grid(D/4, D/4, D/4), &ta, sizeof(ta)));

        tex_read_args_t ra = {0};
        ra.dim[0] = D; ra.dim[1] = D; ra.dim[2] = D;
        ra.tex = tex; ra.dst = DEV_ADDR(d);
        ASSERT_TRUE(md_gpu_launch(f.compute, f.k_tex_read, md_gpu_grid(D/4, D/4, D/4), &ra, sizeof(ra)));

        memset(host, 0, VOX * sizeof(float));
        ASSERT_TRUE(md_gpu_memcpy_async(host, d, VOX * sizeof(float), f.compute));
        md_gpu_stream_sync(f.compute);
        md_gpu_device_poll(f.dev);

        for (int i = 0; i < VOX; i += 61) ASSERT_EQ((float)i * scale, host[i]);
    }

    md_gpu_free(d, f.compute);
    free(host);
    md_gpu_tex_destroy(tex, f.compute);
    md_gpu_device_poll(f.dev);
    gpu_close(&f);
}

/* Many transitions in one command buffer: compute, transfer, compute,
   transfer... Each step depends on the one before, so a single missed
   transition anywhere in the chain shows up at the end. */
UTEST(gpu, hazard_alternating_encoder_chain) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 1 << 18, STEPS = 12 };
    uint32_t* host = (uint32_t*)malloc(N * sizeof(uint32_t));
    ASSERT_TRUE(host != NULL);
    uint32_t* a = (uint32_t*)md_gpu_malloc(f.pool, N * sizeof(uint32_t), f.compute);
    uint32_t* b = (uint32_t*)md_gpu_malloc(f.pool, N * sizeof(uint32_t), f.compute);
    ASSERT_TRUE(a && b);

    fill_args_t fa = {0};
    fa.n = N; fa.base = 0; fa.dst = DEV_ADDR(a);
    ASSERT_TRUE(md_gpu_launch(f.compute, f.k_fill, md_gpu_grid_1d(N, 64), &fa, sizeof(fa)));

    for (int step = 0; step < STEPS; ++step) {
        /* transfer: a -> b */
        ASSERT_TRUE(md_gpu_memcpy_async(b, a, N * sizeof(uint32_t), f.compute));
        /* compute: b -> a, +1 */
        scale_args_t sa = {0};
        sa.n = N; sa.mul = 1; sa.add = 1;
        sa.src = DEV_ADDR(b); sa.dst = DEV_ADDR(a);
        ASSERT_TRUE(md_gpu_launch(f.compute, f.k_scale, md_gpu_grid_1d(N, 64), &sa, sizeof(sa)));
    }

    memset(host, 0, N * sizeof(uint32_t));
    ASSERT_TRUE(md_gpu_memcpy_async(host, a, N * sizeof(uint32_t), f.compute));
    md_gpu_stream_sync(f.compute);
    md_gpu_device_poll(f.dev);

    for (int i = 0; i < N; i += 97) ASSERT_EQ((uint32_t)(i + STEPS), host[i]);

    md_gpu_free(a, f.compute);
    md_gpu_free(b, f.compute);
    free(host);
    gpu_close(&f);
}

/* The same producer/consumer pair, forced into separate submissions. Different
   mechanism again -- a timeline wait rather than an encoder-level one. */
UTEST(gpu, hazard_across_submission_boundary) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    float* host = (float*)malloc(HAZ_VOX * sizeof(float));
    ASSERT_TRUE(host != NULL);

    md_gpu_tex_desc_t td = {0};
    td.width = HAZ_D; td.height = HAZ_D; td.depth = HAZ_D;
    td.format = MD_GPU_FORMAT_R32_FLOAT;
    td.flags  = MD_GPU_TEX_STORAGE;
    md_gpu_tex_t tex = md_gpu_tex_create(f.dev, &td);
    ASSERT_TRUE(tex != 0);

    for (int r = 0; r < HAZ_ROUNDS; ++r) {
        const float scale = (float)(r + 1) * 3.0f;
        tex_args_t ta = {0};
        ta.dim[0] = HAZ_D; ta.dim[1] = HAZ_D; ta.dim[2] = HAZ_D;
        ta.tex = tex; ta.scale = scale;
        ASSERT_TRUE(md_gpu_launch(f.compute, f.k_tex_write,
                                  md_gpu_grid(HAZ_D/4, HAZ_D/4, HAZ_D/4), &ta, sizeof(ta)));

        md_gpu_stream_flush(f.compute);     /* producer and consumer now split */

        memset(host, 0, HAZ_VOX * sizeof(float));
        ASSERT_TRUE(md_gpu_memcpy_from_tex_async(host, tex, NULL, HAZ_VOX * sizeof(float), f.compute));
        md_gpu_stream_sync(f.compute);
        md_gpu_device_poll(f.dev);

        for (int i = 0; i < HAZ_VOX; i += HAZ_STRIDE) ASSERT_EQ((float)i * scale, host[i]);
    }

    free(host);
    md_gpu_tex_destroy(tex, f.compute);
    md_gpu_device_poll(f.dev);
    gpu_close(&f);
}

/* Producer and consumer on different streams, joined by an explicit sync. This
   is the one case where the ordering is the caller's to establish, so it checks
   that record/wait actually carries a texture dependency. */
UTEST(gpu, hazard_across_streams_via_sync) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_stream_t producer = md_gpu_stream_create(f.dev, MD_GPU_STREAM_COMPUTE, "producer");
    ASSERT_TRUE(producer != NULL);

    float* host = (float*)malloc(HAZ_VOX * sizeof(float));
    ASSERT_TRUE(host != NULL);

    md_gpu_tex_desc_t td = {0};
    td.width = HAZ_D; td.height = HAZ_D; td.depth = HAZ_D;
    td.format = MD_GPU_FORMAT_R32_FLOAT;
    td.flags  = MD_GPU_TEX_STORAGE;
    md_gpu_tex_t tex = md_gpu_tex_create(f.dev, &td);
    ASSERT_TRUE(tex != 0);

    for (int r = 0; r < HAZ_ROUNDS; ++r) {
        const float scale = (float)(r + 1) * 7.0f;
        tex_args_t ta = {0};
        ta.dim[0] = HAZ_D; ta.dim[1] = HAZ_D; ta.dim[2] = HAZ_D;
        ta.tex = tex; ta.scale = scale;
        ASSERT_TRUE(md_gpu_launch(producer, f.k_tex_write,
                                  md_gpu_grid(HAZ_D/4, HAZ_D/4, HAZ_D/4), &ta, sizeof(ta)));

        md_gpu_sync_t done = md_gpu_stream_record(producer);
        ASSERT_TRUE(md_gpu_sync_is_valid(done));
        md_gpu_stream_wait(f.compute, done);

        memset(host, 0, HAZ_VOX * sizeof(float));
        ASSERT_TRUE(md_gpu_memcpy_from_tex_async(host, tex, NULL, HAZ_VOX * sizeof(float), f.compute));
        md_gpu_stream_sync(f.compute);
        md_gpu_device_poll(f.dev);

        for (int i = 0; i < HAZ_VOX; i += HAZ_STRIDE) ASSERT_EQ((float)i * scale, host[i]);
    }

    free(host);
    md_gpu_tex_destroy(tex, f.compute);
    md_gpu_device_poll(f.dev);
    md_gpu_stream_destroy(producer);
    gpu_close(&f);
}

/* The same dependency captured into a graph. Replay re-encodes the operations,
   so whatever ordering the backend inserts during normal recording has to be
   inserted again here -- an easy thing to implement once and forget. */
UTEST(gpu, hazard_inside_a_replayed_graph) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { D = 64, VOX = D * D * D, REPLAYS = 4 };
    md_gpu_tex_desc_t td = {0};
    td.width = D; td.height = D; td.depth = D;
    td.format = MD_GPU_FORMAT_R32_FLOAT;
    td.flags  = MD_GPU_TEX_STORAGE;
    md_gpu_tex_t tex = md_gpu_tex_create(f.dev, &td);
    ASSERT_TRUE(tex != 0);

    float* host = (float*)malloc(VOX * sizeof(float));
    ASSERT_TRUE(host != NULL);
    float* d = (float*)md_gpu_malloc(f.pool, VOX * sizeof(float), f.compute);
    ASSERT_TRUE(d != NULL);

    /* write the volume, then read it straight back into a device buffer */
    ASSERT_TRUE(md_gpu_capture_begin(f.compute, "hazard graph"));
    tex_args_t ta = {0};
    ta.dim[0] = D; ta.dim[1] = D; ta.dim[2] = D;
    ta.tex = tex; ta.scale = 2.0f;
    ASSERT_TRUE(md_gpu_launch(f.compute, f.k_tex_write, md_gpu_grid(D/4, D/4, D/4), &ta, sizeof(ta)));
    ASSERT_TRUE(md_gpu_memcpy_from_tex_async(d, tex, NULL, VOX * sizeof(float), f.compute));
    md_gpu_graph_t g = md_gpu_capture_end(f.compute);
    ASSERT_TRUE(g != NULL);

    for (int r = 0; r < REPLAYS; ++r) {
        ASSERT_TRUE(md_gpu_memset_async(d, 0, VOX * sizeof(float), f.compute));
        ASSERT_TRUE(md_gpu_graph_launch(g, f.compute));

        memset(host, 0, VOX * sizeof(float));
        ASSERT_TRUE(md_gpu_memcpy_async(host, d, VOX * sizeof(float), f.compute));
        md_gpu_stream_sync(f.compute);
        md_gpu_device_poll(f.dev);

        for (int i = 0; i < VOX; i += 61) ASSERT_EQ((float)i * 2.0f, host[i]);
    }

    md_gpu_graph_destroy(g);
    md_gpu_free(d, f.compute);
    free(host);
    md_gpu_tex_destroy(tex, f.compute);
    md_gpu_device_poll(f.dev);
    gpu_close(&f);
}

/* Samplers are part of the public surface and never touched by a test. */
UTEST(gpu, sampler_create_and_destroy) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_sampler_desc_t sd = {0};
    sd.min_filter = MD_GPU_FILTER_LINEAR;
    sd.mag_filter = MD_GPU_FILTER_LINEAR;
    sd.address_u = sd.address_v = sd.address_w = MD_GPU_ADDRESS_CLAMP_TO_EDGE;
    sd.label = "linear clamp";
    md_gpu_sampler_t a = md_gpu_sampler_create(f.dev, &sd);
    EXPECT_TRUE(a != 0);

    md_gpu_sampler_t b = md_gpu_sampler_create(f.dev, NULL);   /* defaults */
    EXPECT_TRUE(b != 0);

    md_gpu_sampler_destroy(a);
    md_gpu_sampler_destroy(b);
    md_gpu_sampler_destroy(0);        /* must be a no-op, not a crash */
    gpu_close(&f);
}

/* A texture asked for both usages hands out two handles on backends where
   storage and sampled descriptors are separate arrays, and the same one where
   they are not. Either way the create handle is what identifies it. */
UTEST(gpu, texture_storage_and_sampled_handles) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_tex_desc_t td = {0};
    td.width = 8; td.height = 8; td.depth = 8;
    td.format = MD_GPU_FORMAT_R32_FLOAT;
    td.flags  = MD_GPU_TEX_STORAGE | MD_GPU_TEX_SAMPLED;
    md_gpu_tex_t tex = md_gpu_tex_create(f.dev, &td);
    ASSERT_TRUE(tex != 0);

    md_gpu_tex_t sampled = md_gpu_tex_sampled(tex);
    EXPECT_TRUE(sampled != 0);

    /* Only the create handle identifies the texture. */
    md_gpu_tex_desc_t back = {0};
    ASSERT_TRUE(md_gpu_tex_desc(tex, &back));
    EXPECT_EQ(8u, back.width);
    EXPECT_EQ((uint32_t)(MD_GPU_TEX_STORAGE | MD_GPU_TEX_SAMPLED), back.flags);

    md_gpu_tex_destroy(tex, f.compute);
    md_gpu_device_poll(f.dev);
    gpu_close(&f);
}

/* =========================================================================
   Synchronisation at scale
   ========================================================================= */

/* Both backends keep pending cross-stream waits in a small fixed array and
   sample timelines into a fixed-size snapshot during poll. Fan more streams in
   than those arrays hold and the joining stream must still see every producer's
   writes. */
UTEST(gpu, wide_cross_stream_fan_in) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { STREAMS = 20, N = 128 };
    md_gpu_stream_t s[STREAMS] = {0};
    uint32_t* buf[STREAMS] = {0};

    for (int i = 0; i < STREAMS; ++i) {
        s[i] = md_gpu_stream_create(f.dev, MD_GPU_STREAM_COMPUTE, "producer");
        ASSERT_TRUE(s[i] != NULL);
        buf[i] = (uint32_t*)md_gpu_malloc(f.pool, N * sizeof(uint32_t), s[i]);
        ASSERT_TRUE(buf[i] != NULL);

        fill_args_t fa = {0};
        fa.n = N; fa.base = (uint32_t)(i * 1000); fa.dst = DEV_ADDR(buf[i]);
        ASSERT_TRUE(md_gpu_launch(s[i], f.k_fill, md_gpu_grid_1d(N, 64), &fa, sizeof(fa)));
    }

    /* One consumer waits on all of them, then reads every buffer. */
    for (int i = 0; i < STREAMS; ++i) {
        md_gpu_sync_t sync = md_gpu_stream_record(s[i]);
        ASSERT_TRUE(md_gpu_sync_is_valid(sync));
        md_gpu_stream_wait(f.compute, sync);
    }

    static uint32_t host[STREAMS][N];
    memset(host, 0, sizeof(host));
    for (int i = 0; i < STREAMS; ++i) {
        ASSERT_TRUE(md_gpu_memcpy_async(host[i], buf[i], N * sizeof(uint32_t), f.compute));
    }
    md_gpu_stream_sync(f.compute);
    md_gpu_device_poll(f.dev);

    for (int i = 0; i < STREAMS; ++i) {
        for (int j = 0; j < N; ++j) ASSERT_EQ((uint32_t)(i * 1000 + j), host[i][j]);
    }

    for (int i = 0; i < STREAMS; ++i) {
        md_gpu_free(buf[i], s[i]);
        md_gpu_stream_destroy(s[i]);
    }
    gpu_close(&f);
}

/* The header promises different streams may be driven concurrently from
   different threads, and that malloc/free are thread-safe. Every thread owns
   its stream and its memory, so any failure here is md_gpu's shared state --
   the live-allocation registry above all -- not the test sharing something. */
typedef struct {
    gpu_fixture_t* f;
    int            index;
    bool           ok;
    char           failure[256];
} gpu_thread_ctx_t;

static void gpu_thread_body(void* user) {
    gpu_thread_ctx_t* c = (gpu_thread_ctx_t*)user;
    enum { N = 512, ITERS = 24 };
    md_gpu_stream_t s = md_gpu_stream_create(c->f->dev, MD_GPU_STREAM_COMPUTE, "worker");
    if (!s) { snprintf(c->failure, sizeof(c->failure), "stream creation: %s", gpu_no_device_reason()); return; }

    bool ok = true;
    for (int it = 0; it < ITERS && ok; ++it) {
        uint32_t* d = (uint32_t*)md_gpu_malloc(c->f->pool, N * sizeof(uint32_t), s);
        if (!d) {
            snprintf(c->failure, sizeof(c->failure), "allocation: %s", gpu_no_device_reason());
            ok = false;
            break;
        }

        fill_args_t fa = {0};
        fa.n = N; fa.base = (uint32_t)(c->index * 100000 + it); fa.dst = DEV_ADDR(d);
        if (!md_gpu_launch(s, c->f->k_fill, md_gpu_grid_1d(N, 64), &fa, sizeof(fa))) {
            snprintf(c->failure, sizeof(c->failure), "launch: %s", gpu_no_device_reason());
            ok = false;
            break;
        }

        uint32_t host[N];
        if (!md_gpu_memcpy_async(host, d, sizeof(host), s)) {
            snprintf(c->failure, sizeof(c->failure), "copy: %s", gpu_no_device_reason());
            ok = false;
            break;
        }
        md_gpu_stream_sync(s);

        for (int i = 0; i < N && ok; ++i) {
            if (host[i] != (uint32_t)(c->index * 100000 + it + i)) {
                snprintf(c->failure, sizeof(c->failure), "readback data mismatch");
                ok = false;
            }
        }
        md_gpu_free(d, s);
    }

    md_gpu_stream_destroy(s);
    c->ok = ok;
}

UTEST(gpu, concurrent_streams_from_threads) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { THREADS = 4 };
    gpu_thread_ctx_t ctx[THREADS];
    md_thread_t* th[THREADS];
    for (int i = 0; i < THREADS; ++i) {
        ctx[i].f = &f; ctx[i].index = i; ctx[i].ok = false; ctx[i].failure[0] = '\0';
        th[i] = md_thread_create(gpu_thread_body, &ctx[i]);
        ASSERT_TRUE(th[i] != NULL);
    }
    for (int i = 0; i < THREADS; ++i) {
        md_thread_join(th[i]);
        EXPECT_TRUE_MSG(ctx[i].ok, ctx[i].failure[0] ? ctx[i].failure : "unknown failure");
    }

    md_gpu_device_poll(f.dev);
    gpu_close(&f);
}

/* =========================================================================
   Lifetime
   ========================================================================= */

/* Destroying a pool releases its device memory. Doing it with work still in
   flight is either legal -- in which case it must wait -- or it is not, and the
   header should say so. Today it says nothing, so pin the behaviour down. */
UTEST(gpu, pool_destroy_with_work_in_flight) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_pool_t pool = md_gpu_pool_create(f.dev, &(md_gpu_pool_desc_t){ .flags = MD_GPU_MEM_DEVICE, .label = "in-flight" });
    ASSERT_TRUE(pool != NULL);

    enum { N = 65536 };
    uint32_t* d = (uint32_t*)md_gpu_malloc(pool, N * sizeof(uint32_t), f.compute);
    ASSERT_TRUE(d != NULL);

    /* Enough launches that the queue is genuinely busy when the pool goes. */
    for (int i = 0; i < 32; ++i) {
        fill_args_t fa = {0};
        fa.n = N; fa.base = (uint32_t)i; fa.dst = DEV_ADDR(d);
        ASSERT_TRUE(md_gpu_launch(f.compute, f.k_fill, md_gpu_grid_1d(N, 64), &fa, sizeof(fa)));
    }
    md_gpu_stream_flush(f.compute);

    md_gpu_pool_destroy(pool);        /* must not free memory the GPU is reading */

    md_gpu_stream_sync(f.compute);
    md_gpu_device_poll(f.dev);
    gpu_close(&f);
}

/* The header says device_destroy "destroys the device and everything created
   from it". Leave a pool, a texture, a sampler, a graph and a stream behind and
   let the leak checkers say whether that is true. */
UTEST(gpu, device_destroy_reclaims_undestroyed_objects) {
    gpu_alloc_stats_t stats = {0};
    md_allocator_i alloc = { (md_allocator_o*)&stats, gpu_test_realloc };

    md_gpu_device_desc_t dd = {0};
    dd.alloc = &alloc;
    md_gpu_device_t dev = md_gpu_device_create(&dd);
    if (!dev) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_stream_t s = md_gpu_stream_default(dev, MD_GPU_STREAM_COMPUTE);
    md_gpu_stream_t extra = md_gpu_stream_create(dev, MD_GPU_STREAM_COMPUTE, "leaked");
    ASSERT_TRUE(extra != NULL);

    md_gpu_pool_t pool = md_gpu_pool_create(dev, &(md_gpu_pool_desc_t){ .flags = MD_GPU_MEM_DEVICE, .label = "leaked" });
    ASSERT_TRUE(pool != NULL);
    ASSERT_TRUE(md_gpu_malloc(pool, 64 * 1024, s) != NULL);

    md_gpu_tex_desc_t td = {0};
    td.width = 4; td.height = 4; td.depth = 4;
    td.format = MD_GPU_FORMAT_R32_FLOAT;
    td.flags  = MD_GPU_TEX_STORAGE;
    ASSERT_TRUE(md_gpu_tex_create(dev, &td) != 0);
    ASSERT_TRUE(md_gpu_sampler_create(dev, NULL) != 0);

    /* Nothing above is destroyed by hand. */
    md_gpu_device_destroy(dev);
    EXPECT_EQ(0u, (unsigned)stats.live_bytes);
}

/* =========================================================================
   Edge cases and error reporting
   ========================================================================= */

UTEST(gpu, null_and_zero_arguments_are_tolerated) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    /* Frees and destroys of nothing. */
    md_gpu_free(NULL, f.compute);
    md_gpu_tex_destroy(0, f.compute);
    md_gpu_kernel_destroy(NULL);
    md_gpu_graph_destroy(NULL);
    md_gpu_pool_trim(NULL, 0);

    /* Zero-sized transfers are no-ops, not failures. */
    uint32_t* d = (uint32_t*)md_gpu_malloc(f.pool, 256, f.compute);
    ASSERT_TRUE(d != NULL);
    uint32_t host = 0;
    EXPECT_TRUE(md_gpu_memcpy_async(d, &host, 0, f.compute));
    EXPECT_TRUE(md_gpu_memset_async(d, 0, 0, f.compute));

    /* An empty grid launches nothing and succeeds. */
    fill_args_t fa = {0};
    fa.n = 1; fa.dst = DEV_ADDR(d);
    EXPECT_TRUE(md_gpu_launch(f.compute, f.k_fill, md_gpu_grid(0, 1, 1), &fa, sizeof(fa)));

    /* Queries on a pointer that is not a live allocation. */
    EXPECT_EQ(0u, (unsigned)md_gpu_ptr_size(NULL));
    EXPECT_TRUE(md_gpu_ptr_base(NULL) == NULL);
    EXPECT_TRUE(md_gpu_host_ptr(NULL) == NULL);

    /* Device-local memory is not host-mappable. */
    EXPECT_TRUE(md_gpu_host_ptr(d) == NULL);
    EXPECT_EQ((uint32_t)MD_GPU_MEM_DEVICE, md_gpu_pool_flags(f.pool));

    md_gpu_stream_sync(f.compute);
    md_gpu_free(d, f.compute);
    gpu_close(&f);
}

/* Host-visible memory is dereferenceable, and pointer arithmetic on the device
   pointer must agree with arithmetic on the host view. */
UTEST(gpu, host_pointer_tracks_device_pointer) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 256 };
    uint32_t* d = (uint32_t*)md_gpu_malloc(f.read_pool, N * sizeof(uint32_t), f.compute);
    ASSERT_TRUE(d != NULL);

    uint32_t* h = (uint32_t*)md_gpu_host_ptr(d);
    ASSERT_TRUE(h != NULL);
    uint32_t* h_mid = (uint32_t*)md_gpu_host_ptr(d + N / 2);
    ASSERT_TRUE(h_mid == h + N / 2);

    fill_args_t fa = {0};
    fa.n = N; fa.base = 9; fa.dst = DEV_ADDR(d);
    ASSERT_TRUE(md_gpu_launch(f.compute, f.k_fill, md_gpu_grid_1d(N, 64), &fa, sizeof(fa)));
    md_gpu_stream_sync(f.compute);

    for (int i = 0; i < N; ++i) EXPECT_EQ((uint32_t)(9 + i), h[i]);

    md_gpu_free(d, f.compute);
    gpu_close(&f);
}

/* A failing call must leave a description behind; that is the only diagnostic
   the API offers. */
UTEST(gpu, last_error_describes_the_failure) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    uint32_t host = 0;
    EXPECT_FALSE(md_gpu_memset_async(&host, 0, sizeof(host), f.compute));
    EXPECT_TRUE(md_gpu_last_error() != NULL);

    /* Two uploads open at once on one stream is documented as illegal. */
    uint32_t* d = (uint32_t*)md_gpu_malloc(f.pool, 256, f.compute);
    ASSERT_TRUE(d != NULL);
    void* p = md_gpu_upload_begin(f.compute, d, 256);
    ASSERT_TRUE(p != NULL);
    EXPECT_TRUE(md_gpu_upload_begin(f.compute, d, 256) == NULL);
    EXPECT_TRUE(md_gpu_upload_end(f.compute));

    md_gpu_stream_sync(f.compute);
    md_gpu_free(d, f.compute);
    gpu_close(&f);
}

UTEST(gpu, device_info_is_sane) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_device_info_t info;
    ASSERT_TRUE(md_gpu_device_info(f.dev, &info));
    EXPECT_GT(info.max_threads_per_group, 0u);
    EXPECT_GT(info.preferred_group_multiple, 0u);
    /* Every kernel in the fixture must fit inside the device limit. */
    md_gpu_kernel_info_t ki;
    ASSERT_TRUE(md_gpu_kernel_info(f.k_tex_write, &ki));
    EXPECT_LE(ki.group_size[0] * ki.group_size[1] * ki.group_size[2], info.max_threads_per_group);
    EXPECT_GT(info.timestamp_period_ns_den, 0u);   /* used as a divisor */

    gpu_close(&f);
}

/* An indirect grid derived from a device-side count of zero must dispatch
   nothing rather than fault -- the natural outcome of a compaction pass that
   found no work. */
UTEST(gpu, indirect_launch_with_zero_count) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 64 };
    uint32_t* data  = (uint32_t*)md_gpu_malloc(f.pool, N * sizeof(uint32_t), f.compute);
    uint32_t* count = (uint32_t*)md_gpu_malloc(f.pool, sizeof(uint32_t), f.compute);
    uint32_t* grid  = (uint32_t*)md_gpu_malloc(f.pool, 3 * sizeof(uint32_t), f.compute);
    ASSERT_TRUE(data && count && grid);

    fill_args_t fa = {0};
    fa.n = N; fa.base = 0; fa.dst = DEV_ADDR(data);
    ASSERT_TRUE(md_gpu_launch(f.compute, f.k_fill, md_gpu_grid_1d(N, 64), &fa, sizeof(fa)));
    ASSERT_TRUE(md_gpu_memset_async(count, 0, sizeof(uint32_t), f.compute));

    const uint32_t local[3] = {64, 1, 1};
    ASSERT_TRUE(md_gpu_make_grid(f.compute, grid, count, local));

    bump_args_t ba = {0};
    ba.n = N; ba.delta = 777; ba.dst = DEV_ADDR(data);
    ASSERT_TRUE(md_gpu_launch_indirect(f.compute, f.k_bump, grid, &ba, sizeof(ba)));

    uint32_t host_grid[3] = {1, 1, 1}, host[N];
    ASSERT_TRUE(md_gpu_memcpy_async(host_grid, grid, sizeof(host_grid), f.compute));
    ASSERT_TRUE(md_gpu_memcpy_async(host, data, sizeof(host), f.compute));
    md_gpu_stream_sync(f.compute);

    EXPECT_EQ(0u, host_grid[0]);
    for (int i = 0; i < N; ++i) EXPECT_EQ((uint32_t)i, host[i]);   /* bump did not run */

    md_gpu_free(data,  f.compute);
    md_gpu_free(count, f.compute);
    md_gpu_free(grid,  f.compute);
    gpu_close(&f);
}

/* Blocks freed on one stream must not be handed to another until the first
   stream has actually passed the free point. */
UTEST(gpu, pool_reuse_across_streams_waits_for_completion) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_pool_t pool = md_gpu_pool_create(f.dev, &(md_gpu_pool_desc_t){ .flags = MD_GPU_MEM_DEVICE, .release_threshold = 16 * 1024 * 1024, .label = "cross-stream" });
    ASSERT_TRUE(pool != NULL);
    md_gpu_stream_t other = md_gpu_stream_create(f.dev, MD_GPU_STREAM_COMPUTE, "other");
    ASSERT_TRUE(other != NULL);

    enum { N = 65536 };
    uint32_t* a = (uint32_t*)md_gpu_malloc(pool, N * sizeof(uint32_t), f.compute);
    ASSERT_TRUE(a != NULL);

    for (int i = 0; i < 16; ++i) {
        fill_args_t fa = {0};
        fa.n = N; fa.base = 4242; fa.dst = DEV_ADDR(a);
        ASSERT_TRUE(md_gpu_launch(f.compute, f.k_fill, md_gpu_grid_1d(N, 64), &fa, sizeof(fa)));
    }
    md_gpu_stream_flush(f.compute);
    md_gpu_free(a, f.compute);

    /* `other` has no ordering with f.compute, so the pool may only hand this
       block over once f.compute has actually passed the free -- otherwise the
       write below races the launches above and the readback sees 4242. */
    uint32_t* b = (uint32_t*)md_gpu_malloc(pool, N * sizeof(uint32_t), other);
    ASSERT_TRUE(b != NULL);

    fill_args_t fb = {0};
    fb.n = N; fb.base = 1; fb.dst = DEV_ADDR(b);
    ASSERT_TRUE(md_gpu_launch(other, f.k_fill, md_gpu_grid_1d(N, 64), &fb, sizeof(fb)));

    static uint32_t host[N];
    ASSERT_TRUE(md_gpu_memcpy_async(host, b, sizeof(host), other));
    md_gpu_stream_sync(other);
    md_gpu_stream_sync(f.compute);
    for (int i = 0; i < N; ++i) ASSERT_EQ((uint32_t)(1 + i), host[i]);

    md_gpu_free(b, other);
    md_gpu_stream_destroy(other);
    md_gpu_pool_destroy(pool);
    gpu_close(&f);
}

/* =========================================================================
   Stream lifetime
   ========================================================================= */

/* A pool block records the stream it was freed on so that a later allocation
   can tell whether reuse is safe. Destroying that stream must not leave the
   record pointing at freed memory -- and must not strand the block either: its
   free point is a timeline value the destroyed stream will now never signal,
   so a block left alone would never be reusable again.

   Needs a caching pool. With the default release_threshold of 0 the block is
   handed straight back to the driver on free, so nothing survives to dangle. */
UTEST(gpu, stream_destroy_releases_blocks_it_freed) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_pool_t pool = md_gpu_pool_create(f.dev, &(md_gpu_pool_desc_t){
        .flags = MD_GPU_MEM_DEVICE, .release_threshold = 16 * 1024 * 1024, .label = "caching" });
    ASSERT_TRUE(pool != NULL);

    enum { N = 1024 };
    md_gpu_stream_t worker = md_gpu_stream_create(f.dev, MD_GPU_STREAM_COMPUTE, "worker");
    ASSERT_TRUE(worker != NULL);

    uint32_t* a = (uint32_t*)md_gpu_malloc(pool, N * sizeof(uint32_t), worker);
    ASSERT_TRUE(a != NULL);
    ASSERT_TRUE(md_gpu_memset_async(a, 0, N * sizeof(uint32_t), worker));
    md_gpu_stream_sync(worker);
    md_gpu_free(a, worker);

    md_gpu_pool_stats_t before = {0};
    md_gpu_pool_stats(pool, &before);

    md_gpu_stream_destroy(worker);

    /* Reading b->free_stream here is what the sanitiser caught. */
    uint32_t* b = (uint32_t*)md_gpu_malloc(pool, N * sizeof(uint32_t), f.compute);
    ASSERT_TRUE(b != NULL);

    /* And it must be the cached block, not a fresh one: if the stale free
       point were left in place the block could never satisfy an allocation. */
    md_gpu_pool_stats_t after = {0};
    md_gpu_pool_stats(pool, &after);
    EXPECT_EQ(before.bytes_reserved, after.bytes_reserved);

    fill_args_t fa = {0};
    fa.n = N; fa.base = 7; fa.dst = DEV_ADDR(b);
    ASSERT_TRUE(md_gpu_launch(f.compute, f.k_fill, md_gpu_grid_1d(N, 64), &fa, sizeof(fa)));
    static uint32_t host[N];
    ASSERT_TRUE(md_gpu_memcpy_async(host, b, sizeof(host), f.compute));
    md_gpu_stream_sync(f.compute);
    for (int i = 0; i < N; ++i) ASSERT_EQ((uint32_t)(7 + i), host[i]);

    md_gpu_free(b, f.compute);
    md_gpu_pool_destroy(pool);
    md_gpu_device_poll(f.dev);
    gpu_close(&f);
}

/* A texture destroyed with work outstanding waits on every stream that had any.
   The wait list used to be a fixed eight entries, so past that the image was
   freed while later streams could still be reading it. */
UTEST(gpu, texture_destroy_waits_on_every_stream) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { STREAMS = 12, D = 8 };
    md_gpu_stream_t s[STREAMS];
    for (int i = 0; i < STREAMS; ++i) {
        s[i] = md_gpu_stream_create(f.dev, MD_GPU_STREAM_COMPUTE, "many");
        ASSERT_TRUE(s[i] != NULL);
    }

    md_gpu_tex_desc_t td = {0};
    td.width = D; td.height = D; td.depth = D;
    td.format = MD_GPU_FORMAT_R32_FLOAT;
    td.flags  = MD_GPU_TEX_STORAGE;
    md_gpu_tex_t tex = md_gpu_tex_create(f.dev, &td);
    ASSERT_TRUE(tex != 0);

    /* Give every stream outstanding work touching the texture. */
    for (int i = 0; i < STREAMS; ++i) {
        tex_args_t ta = {0};
        ta.dim[0] = D; ta.dim[1] = D; ta.dim[2] = D;
        ta.tex = tex; ta.scale = 1.0f;
        ASSERT_TRUE(md_gpu_launch(s[i], f.k_tex_write, md_gpu_grid(D/4, D/4, D/4), &ta, sizeof(ta)));
        md_gpu_stream_flush(s[i]);
    }

    md_gpu_tex_destroy(tex, NULL);
    md_gpu_device_poll(f.dev);

    for (int i = 0; i < STREAMS; ++i) md_gpu_stream_destroy(s[i]);
    md_gpu_device_poll(f.dev);
    gpu_close(&f);
}

/* The copy transfers the whole region whatever `size` says, so a buffer too
   small for it has to be rejected rather than overrun. */
UTEST(gpu, texture_copy_rejects_undersized_buffer) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { D = 8, VOXELS = D * D * D, FULL = VOXELS * sizeof(float) };
    md_gpu_tex_desc_t td = {0};
    td.width = D; td.height = D; td.depth = D;
    td.format = MD_GPU_FORMAT_R32_FLOAT;
    td.flags  = MD_GPU_TEX_STORAGE;
    md_gpu_tex_t tex = md_gpu_tex_create(f.dev, &td);
    ASSERT_TRUE(tex != 0);

    void* rb = md_gpu_malloc(f.read_pool, FULL, f.compute);
    ASSERT_TRUE(rb != NULL);

    /* One float short in each direction. */
    EXPECT_FALSE(md_gpu_memcpy_from_tex_async(rb, tex, NULL, FULL - sizeof(float), f.compute));
    static float src[VOXELS];
    EXPECT_FALSE(md_gpu_memcpy_to_tex_async(tex, NULL, src, FULL - sizeof(float), f.compute));

    /* The exact size still works. */
    EXPECT_TRUE(md_gpu_memcpy_to_tex_async(tex, NULL, src, FULL, f.compute));
    EXPECT_TRUE(md_gpu_memcpy_from_tex_async(rb, tex, NULL, FULL, f.compute));
    md_gpu_stream_sync(f.compute);

    md_gpu_free(rb, f.compute);
    md_gpu_tex_destroy(tex, f.compute);
    md_gpu_device_poll(f.dev);
    gpu_close(&f);
}

#endif /* MD_ENABLE_GPU */
