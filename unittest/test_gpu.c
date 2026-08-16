#include "utest.h"

#include <core/md_gpu.h>
#include <core/md_allocator.h>

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
} gpu_fixture_t;

#define GPU_KERNEL(fix, field, sym, name)                                       \
    do {                                                                          \
        md_gpu_kernel_desc_t kd = {0};                                            \
        kd.code      = sym##_start;                                               \
        kd.code_size = sym##_byte_size;                                           \
        kd.label     = name;                                                      \
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

    GPU_KERNEL(f, k_fill,      md_shader_gpu_test_fill,       "fill");
    GPU_KERNEL(f, k_scale,     md_shader_gpu_test_scale_add,  "scale_add");
    GPU_KERNEL(f, k_sum,       md_shader_gpu_test_sum_reduce, "sum_reduce");
    GPU_KERNEL(f, k_tex_write, md_shader_gpu_test_tex_write,  "tex_write");
    GPU_KERNEL(f, k_tex_read,  md_shader_gpu_test_tex_read,   "tex_read");
    GPU_KERNEL(f, k_bump,      md_shader_gpu_test_bump,       "bump");
    return true;
}

static void gpu_close(gpu_fixture_t* f) {
    if (f->k_fill)      md_gpu_kernel_destroy(f->k_fill);
    if (f->k_scale)     md_gpu_kernel_destroy(f->k_scale);
    if (f->k_sum)       md_gpu_kernel_destroy(f->k_sum);
    if (f->k_tex_write) md_gpu_kernel_destroy(f->k_tex_write);
    if (f->k_tex_read)  md_gpu_kernel_destroy(f->k_tex_read);
    if (f->k_bump)      md_gpu_kernel_destroy(f->k_bump);
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

#define DEV_ADDR(p) ((uint64_t)(uintptr_t)(p))

/* =========================================================================
   Device and streams
   ========================================================================= */

UTEST(gpu, device_create_destroy) {
    md_gpu_device_t dev = md_gpu_device_create(NULL);
    if (!dev) UTEST_SKIP("No GPU device available");
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
    if (!dev) UTEST_SKIP("No GPU device available");

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
    if (!gpu_open(&f)) UTEST_SKIP("No GPU device available");

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
    if (!gpu_open(&f)) UTEST_SKIP("No GPU device available");

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
    if (!gpu_open(&f)) UTEST_SKIP("No GPU device available");

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
    if (!gpu_open(&f)) UTEST_SKIP("No GPU device available");

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
    if (!gpu_open(&f)) UTEST_SKIP("No GPU device available");

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
    if (!gpu_open(&f)) UTEST_SKIP("No GPU device available");

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
    if (!gpu_open(&f)) UTEST_SKIP("No GPU device available");

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
    if (!gpu_open(&f)) UTEST_SKIP("No GPU device available");

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
    if (!gpu_open(&f)) UTEST_SKIP("No GPU device available");

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
    if (!gpu_open(&f)) UTEST_SKIP("No GPU device available");

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
    if (!gpu_open(&f)) UTEST_SKIP("No GPU device available");

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
    if (!gpu_open(&f)) UTEST_SKIP("No GPU device available");

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
    if (!gpu_open(&f)) UTEST_SKIP("No GPU device available");

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

    for (int i = 0; i < N; ++i) EXPECT_EQ((uint32_t)i << STEPS, host[i]);

    md_gpu_free(a, f.compute);
    md_gpu_free(b, f.compute);
    gpu_close(&f);
}

UTEST(gpu, cross_stream_dependency) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP("No GPU device available");

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
    if (!gpu_open(&f)) UTEST_SKIP("No GPU device available");

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
    if (!gpu_open(&f)) UTEST_SKIP("No GPU device available");

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
    if (!gpu_open(&f)) UTEST_SKIP("No GPU device available");

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
    if (!gpu_open(&f)) UTEST_SKIP("No GPU device available");

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
    if (!gpu_open(&f)) UTEST_SKIP("No GPU device available");

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
    if (!gpu_open(&f)) UTEST_SKIP("No GPU device available");

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
    if (!gpu_open(&f)) UTEST_SKIP("No GPU device available");

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
    if (!gpu_open(&f)) UTEST_SKIP("No GPU device available");

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
    if (!gpu_open(&f)) UTEST_SKIP("No GPU device available");

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
    if (!gpu_open(&f)) UTEST_SKIP("No GPU device available");

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
    if (!gpu_open(&f)) UTEST_SKIP("No GPU device available");

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
    if (!gpu_open(&f)) UTEST_SKIP("No GPU device available");

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
    if (!gpu_open(&f)) UTEST_SKIP("No GPU device available");

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

#endif /* MD_ENABLE_GPU */
