#include "utest.h"

#include <core/md_gpu.h>
#include <core/md_allocator.h>
#include <core/md_os.h>

#include <stdint.h>
#include <stdio.h>
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
    md_gpu_pool_t   write_pool; /* host-writable */
    md_gpu_kernel_t k_fill;
    md_gpu_kernel_t k_scale;
    md_gpu_kernel_t k_sum;
    md_gpu_kernel_t k_tex_write;
    md_gpu_kernel_t k_tex_read;
    md_gpu_kernel_t k_bump;
    md_gpu_kernel_t k_layout;
    md_gpu_kernel_t k_tex_probe;
    md_gpu_kernel_t k_sample;
    md_gpu_kernel_t k_spin;
} gpu_fixture_t;

/* Kernels come from the generated descriptors, which carry the group size and
   argument-struct size read out of the compiled shader. Nothing here repeats
   [numthreads]. */
#define GPU_KERNEL(fix, field, entry)                                          \
    do {                                                                       \
        md_gpu_kernel_desc_t kd = md_shader_gpu_test_##entry##_kernel();      \
        (fix)->field = md_gpu_kernel_create((fix)->dev, &kd);                 \
    } while (0)

static md_gpu_pool_t gpu_pool(md_gpu_device_t dev, md_gpu_mem_kind_t kind, const char* label) {
    md_gpu_pool_desc_t pd = {0};
    pd.kind  = kind;
    pd.label = label;
    return md_gpu_pool_create(dev, &pd);
}

static bool gpu_open(gpu_fixture_t* f) {
    memset(f, 0, sizeof(*f));
    md_gpu_device_desc_t dd = {0};
    dd.enable_validation = true;
    dd.label             = "md_gpu unittest";
    f->dev = md_gpu_device_create(&dd);
    if (!f->dev) return false;

    f->compute    = md_gpu_stream_default(f->dev, MD_GPU_STREAM_COMPUTE);
    f->transfer   = md_gpu_stream_default(f->dev, MD_GPU_STREAM_TRANSFER);
    f->pool       = gpu_pool(f->dev, MD_GPU_MEM_DEVICE,     "test device");
    f->read_pool  = gpu_pool(f->dev, MD_GPU_MEM_HOST_READ,  "test readback");
    f->write_pool = gpu_pool(f->dev, MD_GPU_MEM_HOST_WRITE, "test upload");
    if (!f->pool || !f->read_pool || !f->write_pool) return false;

    GPU_KERNEL(f, k_fill,      fill);
    GPU_KERNEL(f, k_scale,     scale_add);
    GPU_KERNEL(f, k_sum,       sum_reduce);
    GPU_KERNEL(f, k_tex_write, tex_write);
    GPU_KERNEL(f, k_tex_read,  tex_read);
    GPU_KERNEL(f, k_bump,      bump);
    GPU_KERNEL(f, k_layout,    layout_probe);
    GPU_KERNEL(f, k_tex_probe, tex_probe);
    GPU_KERNEL(f, k_sample,    sample_read);
    GPU_KERNEL(f, k_spin,      spin);
    return f->k_fill && f->k_scale && f->k_sum && f->k_tex_write && f->k_tex_read &&
           f->k_bump && f->k_layout && f->k_tex_probe && f->k_sample && f->k_spin;
}

/* Why the fixture could not be opened -- no Vulkan loader, no driver (ICD), no
   compute queue, a kernel that failed to build. Without this a CI log shows
   only "skipped" and gives no way to tell a missing GPU from a broken build. */
static const char* gpu_no_device_reason(void) {
    const char* err = md_gpu_last_error();
    return (err && err[0]) ? err : "No GPU device available";
}

static void gpu_close(gpu_fixture_t* f) {
    md_gpu_kernel_t* ks[] = { &f->k_fill, &f->k_scale, &f->k_sum, &f->k_tex_write, &f->k_tex_read,
                              &f->k_bump, &f->k_layout, &f->k_tex_probe, &f->k_sample, &f->k_spin };
    for (size_t i = 0; i < sizeof(ks) / sizeof(ks[0]); ++i) md_gpu_kernel_destroy(*ks[i]);
    md_gpu_pool_destroy(f->pool);
    md_gpu_pool_destroy(f->read_pool);
    md_gpu_pool_destroy(f->write_pool);
    md_gpu_device_destroy(f->dev);
}

/* Device-local allocation from the fixture pool. */
static md_gpu_addr_t gpu_alloc(gpu_fixture_t* f, md_gpu_stream_t s, size_t size) {
    return md_gpu_malloc(s, f->pool, size).gpu;
}

/* The one way results reach the host: copy into HOST_READ memory, wait, read
   its CPU pointer. Synchronises `s`. */
static bool gpu_read(gpu_fixture_t* f, md_gpu_stream_t s, void* dst, md_gpu_addr_t src, size_t size) {
    md_gpu_mem_t rb = md_gpu_malloc(s, f->read_pool, size);
    if (!rb.cpu) return false;
    bool ok = md_gpu_copy(s, rb.gpu, src, size);
    md_gpu_stream_sync(s);
    if (ok) memcpy(dst, rb.cpu, size);
    md_gpu_free(s, rb.gpu);
    return ok;
}

static bool gpu_read_tex(gpu_fixture_t* f, md_gpu_stream_t s, void* dst, md_gpu_texture_t t, const md_gpu_tex_region_t* region) {
    size_t size = md_gpu_texture_region_size(t, region);
    if (size == 0) return false;
    md_gpu_mem_t rb = md_gpu_malloc(s, f->read_pool, size);
    if (!rb.cpu) return false;
    bool ok = md_gpu_copy_from_texture(s, rb.gpu, t, region);
    md_gpu_stream_sync(s);
    if (ok) memcpy(dst, rb.cpu, size);
    md_gpu_free(s, rb.gpu);
    return ok;
}

/* A cubic R32_FLOAT 3D texture from the fixture pool, created on `s`. */
static md_gpu_texture_t gpu_volume(gpu_fixture_t* f, md_gpu_stream_t s, uint32_t d, md_gpu_tex_usage_t usage) {
    md_gpu_texture_desc_t td = {0};
    td.type   = MD_GPU_TEX_3D;
    td.format = MD_GPU_FORMAT_R32_FLOAT;
    td.usage  = usage;
    td.width  = d; td.height = d; td.depth_or_layers = d;
    td.label  = "volume";
    return md_gpu_texture_create(s, f->pool, &td);
}

/* Argument structs, mirroring unittest/shaders/gpu_test.slang. Pointer fields
   are md_gpu_addr_t, so they take an allocation's .gpu with no cast. */
typedef struct { uint32_t n, base, pad0, pad1; md_gpu_addr_t dst; }               fill_args_t;
typedef struct { uint32_t n, mul, add, pad; md_gpu_addr_t src, dst; }             scale_args_t;
typedef struct { uint32_t n, pad0, pad1, pad2; md_gpu_addr_t src, out_val; }      sum_args_t;
typedef struct { uint32_t dim[3]; float scale; md_gpu_storage_tex_t tex; }        tex_args_t;
typedef struct { uint32_t dim[3], pad; md_gpu_storage_tex_t tex; md_gpu_addr_t dst; } tex_read_args_t;
typedef struct { uint32_t n, delta, pad0, pad1; md_gpu_addr_t dst; }              bump_args_t;
typedef struct { uint32_t dim[3], pad; md_gpu_sampled_tex_t tex; md_gpu_sampler_t smp; md_gpu_addr_t dst; } sample_args_t;
typedef struct { uint32_t n, iters, pad0, pad1; md_gpu_addr_t dst; }              spin_args_t;

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
    md_gpu_addr_t   dst;
} layout_probe_args_t;

typedef struct {
    uint32_t             n, pad0, pad1, pad2;
    md_gpu_addr_t        dst;
    md_gpu_storage_tex_t tex;
    uint32_t             marker, pad3;
} tex_probe_args_t;

static md_gpu_grid_t grid1(md_gpu_kernel_t k, uint32_t n) { return md_gpu_grid_for(k, n, 1, 1); }

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
    md_gpu_pool_t pool = gpu_pool(dev, MD_GPU_MEM_DEVICE, "routes");
    ASSERT_TRUE(pool != NULL);
    md_gpu_mem_t m = md_gpu_malloc(s, pool, 4096);
    ASSERT_TRUE(m.gpu != 0);
    md_gpu_free(s, m.gpu);
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

    /* A stream that has never submitted reports a none sync and syncs instantly. */
    md_gpu_sync_t none = md_gpu_stream_record(a);
    EXPECT_FALSE(md_gpu_sync_is_valid(none));
    EXPECT_TRUE(md_gpu_sync_is_complete(none));
    md_gpu_stream_sync(a);

    md_gpu_stream_destroy(a);
    gpu_close(&f);
}

/* With nothing new issued, record returns the previous submission -- still a
   correct "everything so far" point -- rather than none. */
UTEST(gpu, record_with_nothing_pending_returns_last_submission) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_addr_t d = gpu_alloc(&f, f.compute, 256);
    ASSERT_TRUE(d != 0);
    ASSERT_TRUE(md_gpu_memset(f.compute, d, 0, 256));
    md_gpu_sync_t first = md_gpu_stream_record(f.compute);
    ASSERT_TRUE(md_gpu_sync_is_valid(first));
    md_gpu_sync_t again = md_gpu_stream_record(f.compute);
    EXPECT_TRUE(again.stream == first.stream);
    EXPECT_EQ(first.value, again.value);

    md_gpu_stream_sync(f.compute);
    md_gpu_free(f.compute, d);
    gpu_close(&f);
}

/* =========================================================================
   Memory
   ========================================================================= */

UTEST(gpu, upload_and_read_roundtrip) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 1024 };
    uint32_t src[N], dst[N];
    for (int i = 0; i < N; ++i) src[i] = (uint32_t)i * 2654435761u;
    memset(dst, 0, sizeof(dst));

    md_gpu_addr_t d = gpu_alloc(&f, f.compute, sizeof(src));
    ASSERT_TRUE(d != 0);
    ASSERT_TRUE(md_gpu_upload(f.compute, d, src, sizeof(src)));
    ASSERT_TRUE(gpu_read(&f, f.compute, dst, d, sizeof(dst)));

    for (int i = 0; i < N; ++i) EXPECT_EQ(src[i], dst[i]);

    md_gpu_free(f.compute, d);
    gpu_close(&f);
}

/* The upload's source is consumed before md_gpu_upload returns. */
UTEST(gpu, upload_source_may_be_reused_immediately) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 256 };
    uint32_t src[N];
    md_gpu_addr_t a = gpu_alloc(&f, f.compute, sizeof(src));
    md_gpu_addr_t b = gpu_alloc(&f, f.compute, sizeof(src));
    ASSERT_TRUE(a && b);
    for (int i = 0; i < N; ++i) src[i] = 111;
    ASSERT_TRUE(md_gpu_upload(f.compute, a, src, sizeof(src)));
    for (int i = 0; i < N; ++i) src[i] = 222;          /* overwrite before any sync */
    ASSERT_TRUE(md_gpu_upload(f.compute, b, src, sizeof(src)));

    uint32_t ha[N], hb[N];
    ASSERT_TRUE(gpu_read(&f, f.compute, ha, a, sizeof(ha)));
    ASSERT_TRUE(gpu_read(&f, f.compute, hb, b, sizeof(hb)));
    for (int i = 0; i < N; ++i) { EXPECT_EQ(111u, ha[i]); EXPECT_EQ(222u, hb[i]); }

    md_gpu_free(f.compute, a);
    md_gpu_free(f.compute, b);
    gpu_close(&f);
}

UTEST(gpu, copy_device_to_device) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 256 };
    uint32_t src[N], dst[N];
    for (int i = 0; i < N; ++i) src[i] = (uint32_t)(i + 7) * 11u;
    memset(dst, 0, sizeof(dst));

    md_gpu_addr_t a = gpu_alloc(&f, f.compute, sizeof(src));
    md_gpu_addr_t b = gpu_alloc(&f, f.compute, sizeof(src));
    ASSERT_TRUE(a && b);

    ASSERT_TRUE(md_gpu_upload(f.compute, a, src, sizeof(src)));
    ASSERT_TRUE(md_gpu_copy(f.compute, b, a, sizeof(src)));
    ASSERT_TRUE(gpu_read(&f, f.compute, dst, b, sizeof(dst)));

    for (int i = 0; i < N; ++i) EXPECT_EQ(src[i], dst[i]);

    md_gpu_free(f.compute, a);
    md_gpu_free(f.compute, b);
    gpu_close(&f);
}

UTEST(gpu, memset_aligned_and_unaligned) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { BYTES = 256 };
    uint8_t host[BYTES];
    md_gpu_addr_t d = gpu_alloc(&f, f.compute, BYTES);
    ASSERT_TRUE(d != 0);

    /* Whole buffer, 4-byte aligned. */
    ASSERT_TRUE(md_gpu_memset(f.compute, d, 0xAB, BYTES));
    memset(host, 0, BYTES);
    ASSERT_TRUE(gpu_read(&f, f.compute, host, d, BYTES));
    for (int i = 0; i < BYTES; ++i) EXPECT_EQ(0xAB, host[i]);

    /* Unaligned offset and length, exercising the staged head/tail path. */
    ASSERT_TRUE(md_gpu_memset(f.compute, d + 3, 0x5C, 10));
    /* And one entirely inside a single word. */
    ASSERT_TRUE(md_gpu_memset(f.compute, d + 101, 0x11, 2));
    memset(host, 0, BYTES);
    ASSERT_TRUE(gpu_read(&f, f.compute, host, d, BYTES));
    for (int i = 0; i < BYTES; ++i) {
        uint8_t expect = (i >= 3 && i < 13) ? 0x5C : (i >= 101 && i < 103) ? 0x11 : 0xAB;
        EXPECT_EQ(expect, host[i]);
    }

    md_gpu_free(f.compute, d);
    gpu_close(&f);
}

/* Addresses are plain byte arithmetic: base + offset is a valid sub-range. */
UTEST(gpu, address_arithmetic_subranges) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 128 };
    md_gpu_addr_t base = gpu_alloc(&f, f.compute, N * 2 * sizeof(uint32_t));
    ASSERT_TRUE(base != 0);
    md_gpu_addr_t upper = base + N * sizeof(uint32_t);

    fill_args_t lo = {0}; lo.n = N; lo.base = 1000; lo.dst = base;
    fill_args_t hi = {0}; hi.n = N; hi.base = 5000; hi.dst = upper;
    ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_fill, grid1(f.k_fill, N), lo));
    ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_fill, grid1(f.k_fill, N), hi));

    uint32_t host[N * 2];
    ASSERT_TRUE(gpu_read(&f, f.compute, host, base, sizeof(host)));
    for (int i = 0; i < N; ++i) {
        EXPECT_EQ((uint32_t)(1000 + i), host[i]);
        EXPECT_EQ((uint32_t)(5000 + i), host[N + i]);
    }

    /* A sub-range can be copied on its own, too. */
    uint32_t half[N];
    ASSERT_TRUE(gpu_read(&f, f.compute, half, upper, sizeof(half)));
    for (int i = 0; i < N; ++i) EXPECT_EQ((uint32_t)(5000 + i), half[i]);

    md_gpu_free(f.compute, base);
    gpu_close(&f);
}

/* Host-visible allocations hand back a CPU pointer; device-local ones do not. */
UTEST(gpu, host_visible_memory_has_a_cpu_pointer) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 256 };
    md_gpu_mem_t dev_mem = md_gpu_malloc(f.compute, f.pool, N * sizeof(uint32_t));
    md_gpu_mem_t rd_mem  = md_gpu_malloc(f.compute, f.read_pool, N * sizeof(uint32_t));
    md_gpu_mem_t wr_mem  = md_gpu_malloc(f.compute, f.write_pool, N * sizeof(uint32_t));
    ASSERT_TRUE(dev_mem.gpu && rd_mem.gpu && wr_mem.gpu);
    EXPECT_TRUE(dev_mem.cpu == NULL);
    ASSERT_TRUE(rd_mem.cpu != NULL);
    ASSERT_TRUE(wr_mem.cpu != NULL);

    /* CPU writes into HOST_WRITE memory are visible to a kernel. */
    uint32_t* w = (uint32_t*)wr_mem.cpu;
    for (int i = 0; i < N; ++i) w[i] = (uint32_t)(i * 5);
    scale_args_t sa = {0};
    sa.n = N; sa.mul = 1; sa.add = 3; sa.src = wr_mem.gpu; sa.dst = rd_mem.gpu;
    ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_scale, grid1(f.k_scale, N), sa));
    md_gpu_stream_sync(f.compute);

    /* And the kernel's writes into HOST_READ memory are visible to the CPU. */
    const uint32_t* r = (const uint32_t*)rd_mem.cpu;
    for (int i = 0; i < N; ++i) EXPECT_EQ((uint32_t)(i * 5 + 3), r[i]);

    EXPECT_EQ((int)MD_GPU_MEM_HOST_READ, (int)md_gpu_pool_kind(f.read_pool));

    md_gpu_free(f.compute, dev_mem.gpu);
    md_gpu_free(f.compute, rd_mem.gpu);
    md_gpu_free(f.compute, wr_mem.gpu);
    gpu_close(&f);
}

UTEST(gpu, upload_begin_end_zero_copy) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 512 };
    md_gpu_addr_t d = gpu_alloc(&f, f.compute, N * sizeof(uint32_t));
    ASSERT_TRUE(d != 0);

    uint32_t* p = (uint32_t*)md_gpu_upload_begin(f.compute, d, N * sizeof(uint32_t));
    ASSERT_TRUE(p != NULL);
    for (int i = 0; i < N; ++i) p[i] = (uint32_t)(i * 3 + 1);
    ASSERT_TRUE(md_gpu_upload_end(f.compute));

    uint32_t host[N];
    ASSERT_TRUE(gpu_read(&f, f.compute, host, d, sizeof(host)));
    for (int i = 0; i < N; ++i) EXPECT_EQ((uint32_t)(i * 3 + 1), host[i]);

    md_gpu_free(f.compute, d);
    gpu_close(&f);
}

/* Bigger than the backends' transient page size, so upload_begin has to commit
   a fresh page rather than carve one up. */
UTEST(gpu, upload_larger_than_one_arena_page) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 512 * 1024 };          /* 2 MiB, well past a 256 KiB page */
    md_gpu_addr_t d = gpu_alloc(&f, f.compute, N * sizeof(uint32_t));
    ASSERT_TRUE(d != 0);

    uint32_t* p = (uint32_t*)md_gpu_upload_begin(f.compute, d, N * sizeof(uint32_t));
    ASSERT_TRUE(p != NULL);
    for (int i = 0; i < N; ++i) p[i] = (uint32_t)i * 2654435761u;
    ASSERT_TRUE(md_gpu_upload_end(f.compute));

    uint32_t* host = (uint32_t*)malloc(N * sizeof(uint32_t));
    ASSERT_TRUE(host != NULL);
    ASSERT_TRUE(gpu_read(&f, f.compute, host, d, N * sizeof(uint32_t)));
    for (int i = 0; i < N; ++i) ASSERT_EQ((uint32_t)i * 2654435761u, host[i]);

    free(host);
    md_gpu_free(f.compute, d);
    gpu_close(&f);
}

/* =========================================================================
   Pools
   ========================================================================= */

UTEST(gpu, pool_reuses_freed_blocks) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_pool_t pool = gpu_pool(f.dev, MD_GPU_MEM_DEVICE, "reuse");
    ASSERT_TRUE(pool != NULL);

    md_gpu_pool_stats_t st = {0};
    md_gpu_addr_t a = md_gpu_malloc(f.compute, pool, 4096).gpu;
    ASSERT_TRUE(a != 0);
    md_gpu_pool_stats(pool, &st);
    EXPECT_EQ(4096u, (unsigned)st.bytes_in_use);
    EXPECT_EQ(1u, st.blocks_in_use);
    EXPECT_EQ(1u, (unsigned)st.alloc_count);
    EXPECT_EQ(0u, (unsigned)st.reuse_count);
    uint64_t reserved = st.bytes_reserved;

    md_gpu_free(f.compute, a);
    md_gpu_pool_stats(pool, &st);
    EXPECT_EQ(0u, (unsigned)st.bytes_in_use);
    EXPECT_EQ(1u, st.blocks_cached);        /* cache_limit 0: keeps everything */

    /* Freed on this stream, so reuse is immediate and hits the same block. */
    md_gpu_addr_t b = md_gpu_malloc(f.compute, pool, 4096).gpu;
    EXPECT_TRUE(b == a);

    md_gpu_pool_stats(pool, &st);
    EXPECT_EQ(reserved, st.bytes_reserved);   /* no new memory was committed */
    EXPECT_EQ(1u, (unsigned)st.reuse_count);
    EXPECT_EQ(4096u, (unsigned)st.bytes_peak_in_use);

    md_gpu_free(f.compute, b);
    md_gpu_pool_destroy(pool);
    gpu_close(&f);
}

/* cache_limit bounds what a pool keeps after md_gpu_free. */
UTEST(gpu, pool_cache_limit_releases_excess) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_pool_desc_t pd = {0};
    pd.kind = MD_GPU_MEM_DEVICE; pd.cache_limit = 8192; pd.label = "limited";
    md_gpu_pool_t pool = md_gpu_pool_create(f.dev, &pd);
    ASSERT_TRUE(pool != NULL);

    enum { N = 8 };
    md_gpu_addr_t a[N];
    for (int i = 0; i < N; ++i) { a[i] = md_gpu_malloc(f.compute, pool, 4096).gpu; ASSERT_TRUE(a[i] != 0); }
    md_gpu_stream_sync(f.compute);
    for (int i = 0; i < N; ++i) md_gpu_free(f.compute, a[i]);

    md_gpu_pool_stats_t st = {0};
    md_gpu_pool_stats(pool, &st);
    EXPECT_LE(st.bytes_cached, (uint64_t)8192);

    md_gpu_pool_destroy(pool);
    gpu_close(&f);
}

/* The CPU arena-reset pattern: drop every allocation in one call, keep the
   memory, and let the next round be served entirely from cache. */
UTEST(gpu, pool_reset_recycles_everything) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_pool_t pool = gpu_pool(f.dev, MD_GPU_MEM_DEVICE, "reset");
    ASSERT_TRUE(pool != NULL);

    enum { N = 8 };
    for (int i = 0; i < N; ++i) ASSERT_TRUE(md_gpu_malloc(f.compute, pool, 4096).gpu != 0);
    md_gpu_texture_t tex = NULL;
    {
        md_gpu_texture_desc_t td = {0};
        td.type = MD_GPU_TEX_2D; td.format = MD_GPU_FORMAT_R32_FLOAT; td.usage = MD_GPU_TEX_STORAGE;
        td.width = 16; td.height = 16;
        tex = md_gpu_texture_create(f.compute, pool, &td);
        ASSERT_TRUE(tex != NULL);
    }

    md_gpu_pool_stats_t st = {0};
    md_gpu_pool_stats(pool, &st);
    EXPECT_GE(st.bytes_in_use, (uint64_t)(N * 4096));
    EXPECT_EQ((uint32_t)N, st.blocks_in_use);

    /* One call frees them all -- the texture included -- keeping the blocks. */
    md_gpu_pool_reset(f.compute, pool);

    md_gpu_pool_stats(pool, &st);
    EXPECT_EQ(0u, (unsigned)st.bytes_in_use);
    EXPECT_EQ(0u, st.blocks_in_use);
    EXPECT_EQ((uint32_t)N, st.blocks_cached);
    uint64_t reserved_before = st.bytes_reserved;

    /* The next round is served entirely from cache and commits nothing new. */
    uint64_t reuse_before = st.reuse_count;
    for (int i = 0; i < N; ++i) ASSERT_TRUE(md_gpu_malloc(f.compute, pool, 4096).gpu != 0);
    md_gpu_pool_stats(pool, &st);
    EXPECT_EQ(reserved_before, st.bytes_reserved);
    EXPECT_EQ(reuse_before + N, st.reuse_count);

    /* Reset again, then hand the memory back for real. */
    md_gpu_pool_reset(f.compute, pool);
    md_gpu_stream_sync(f.compute);
    md_gpu_pool_trim(pool, 0);
    md_gpu_pool_stats(pool, &st);
    EXPECT_EQ(0u, (unsigned)st.bytes_reserved);

    md_gpu_pool_destroy(pool);
    md_gpu_device_poll(f.dev);
    gpu_close(&f);
}

/* Reset while the GPU is still using the memory must not disturb it: the
   blocks only become reusable at that point in the stream. */
UTEST(gpu, pool_reset_is_stream_ordered) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_pool_t pool = gpu_pool(f.dev, MD_GPU_MEM_DEVICE, "reset-order");
    ASSERT_TRUE(pool != NULL);

    enum { N = 1024 };
    md_gpu_addr_t d = md_gpu_malloc(f.compute, pool, N * sizeof(uint32_t)).gpu;
    ASSERT_TRUE(d != 0);
    md_gpu_mem_t rb = md_gpu_malloc(f.compute, f.read_pool, N * sizeof(uint32_t));
    ASSERT_TRUE(rb.cpu != NULL);

    fill_args_t a = {0};
    a.n = N; a.base = 5; a.dst = d;
    ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_fill, grid1(f.k_fill, N), a));
    ASSERT_TRUE(md_gpu_copy(f.compute, rb.gpu, d, N * sizeof(uint32_t)));

    /* Reset with the launch and the readback still in flight. */
    md_gpu_pool_reset(f.compute, pool);

    md_gpu_stream_sync(f.compute);
    const uint32_t* host = (const uint32_t*)rb.cpu;
    for (int i = 0; i < N; ++i) EXPECT_EQ((uint32_t)(5 + i), host[i]);

    md_gpu_free(f.compute, rb.gpu);
    md_gpu_pool_destroy(pool);
    gpu_close(&f);
}

UTEST(gpu, pool_trim_releases) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_pool_t pool = gpu_pool(f.dev, MD_GPU_MEM_DEVICE, "trim");
    ASSERT_TRUE(pool != NULL);

    md_gpu_addr_t p = md_gpu_malloc(f.compute, pool, 1024 * 1024).gpu;
    ASSERT_TRUE(p != 0);
    md_gpu_free(f.compute, p);
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

/* Blocks freed on one stream must not be handed to another until the first
   stream has actually passed the free point -- and malloc must never wait for
   that: it commits a new block instead. */
UTEST(gpu, pool_reuse_across_streams_waits_for_completion) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_pool_t pool = gpu_pool(f.dev, MD_GPU_MEM_DEVICE, "cross-stream");
    ASSERT_TRUE(pool != NULL);
    md_gpu_stream_t other = md_gpu_stream_create(f.dev, MD_GPU_STREAM_COMPUTE, "other");
    ASSERT_TRUE(other != NULL);

    enum { N = 65536 };
    md_gpu_addr_t a = md_gpu_malloc(f.compute, pool, N * sizeof(uint32_t)).gpu;
    ASSERT_TRUE(a != 0);

    for (int i = 0; i < 16; ++i) {
        fill_args_t fa = {0};
        fa.n = N; fa.base = 4242; fa.dst = a;
        ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_fill, grid1(f.k_fill, N), fa));
    }
    md_gpu_stream_flush(f.compute);
    md_gpu_free(f.compute, a);

    md_gpu_addr_t b = md_gpu_malloc(other, pool, N * sizeof(uint32_t)).gpu;
    ASSERT_TRUE(b != 0);

    fill_args_t fb = {0};
    fb.n = N; fb.base = 1; fb.dst = b;
    ASSERT_TRUE(MD_GPU_LAUNCH(other, f.k_fill, grid1(f.k_fill, N), fb));

    uint32_t* host = (uint32_t*)malloc(N * sizeof(uint32_t));
    ASSERT_TRUE(host != NULL);
    ASSERT_TRUE(gpu_read(&f, other, host, b, N * sizeof(uint32_t)));
    md_gpu_stream_sync(f.compute);
    for (int i = 0; i < N; ++i) ASSERT_EQ((uint32_t)(1 + i), host[i]);
    free(host);

    md_gpu_free(other, b);
    md_gpu_stream_destroy(other);
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
    md_gpu_addr_t d = gpu_alloc(&f, f.compute, N * sizeof(uint32_t));
    ASSERT_TRUE(d != 0);

    fill_args_t a = {0};
    a.n = N; a.base = 100; a.dst = d;
    ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_fill, grid1(f.k_fill, N), a));

    uint32_t host[N];
    ASSERT_TRUE(gpu_read(&f, f.compute, host, d, sizeof(host)));
    for (int i = 0; i < N; ++i) EXPECT_EQ((uint32_t)(100 + i), host[i]);

    md_gpu_free(f.compute, d);
    gpu_close(&f);
}

/* The central guarantee: consecutive launches in one stream see each other's
   writes with no barrier, no usage declaration and no fence from the caller. */
UTEST(gpu, program_order_chain) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 4096, STEPS = 40 };
    md_gpu_addr_t a = gpu_alloc(&f, f.compute, N * sizeof(uint32_t));
    md_gpu_addr_t b = gpu_alloc(&f, f.compute, N * sizeof(uint32_t));
    ASSERT_TRUE(a && b);

    fill_args_t fa = {0};
    fa.n = N; fa.base = 0; fa.dst = a;
    ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_fill, grid1(f.k_fill, N), fa));

    md_gpu_addr_t src = a, dst = b;
    for (int step = 0; step < STEPS; ++step) {
        scale_args_t sa = {0};
        sa.n = N; sa.mul = 1; sa.add = 1; sa.src = src; sa.dst = dst;
        ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_scale, grid1(f.k_scale, N), sa));
        md_gpu_addr_t tmp = src; src = dst; dst = tmp;
    }

    uint32_t host[N];
    ASSERT_TRUE(gpu_read(&f, f.compute, host, src, sizeof(host)));
    for (int i = 0; i < N; ++i) EXPECT_EQ((uint32_t)(i + STEPS), host[i]);

    md_gpu_free(f.compute, a);
    md_gpu_free(f.compute, b);
    gpu_close(&f);
}

/* The same chain, split across several submissions: the dependency must
   survive command-buffer and submission boundaries. */
UTEST(gpu, program_order_across_submissions) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 1024, STEPS = 8 };
    md_gpu_addr_t a = gpu_alloc(&f, f.compute, N * sizeof(uint32_t));
    md_gpu_addr_t b = gpu_alloc(&f, f.compute, N * sizeof(uint32_t));
    ASSERT_TRUE(a && b);

    fill_args_t fa = {0};
    fa.n = N; fa.base = 0; fa.dst = a;
    ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_fill, grid1(f.k_fill, N), fa));

    md_gpu_addr_t src = a, dst = b;
    for (int step = 0; step < STEPS; ++step) {
        scale_args_t sa = {0};
        sa.n = N; sa.mul = 2; sa.add = 0; sa.src = src; sa.dst = dst;
        ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_scale, grid1(f.k_scale, N), sa));
        md_gpu_stream_flush(f.compute);          /* force a submission boundary */
        md_gpu_addr_t tmp = src; src = dst; dst = tmp;
    }

    uint32_t host[N];
    ASSERT_TRUE(gpu_read(&f, f.compute, host, src, sizeof(host)));
    for (int i = 0; i < N; ++i) EXPECT_EQ(((uint32_t)i << STEPS), host[i]);

    md_gpu_free(f.compute, a);
    md_gpu_free(f.compute, b);
    gpu_close(&f);
}

UTEST(gpu, cross_stream_dependency) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 2048 };
    uint32_t src[N];
    for (int i = 0; i < N; ++i) src[i] = (uint32_t)i;

    md_gpu_addr_t d   = gpu_alloc(&f, f.transfer, sizeof(src));
    md_gpu_addr_t out = gpu_alloc(&f, f.transfer, sizeof(src));
    ASSERT_TRUE(d && out);

    /* Upload on the transfer stream. */
    ASSERT_TRUE(md_gpu_upload(f.transfer, d, src, sizeof(src)));
    md_gpu_sync_t uploaded = md_gpu_stream_record(f.transfer);

    /* Compute stream waits for it, then consumes the data. */
    md_gpu_stream_wait(f.compute, uploaded);
    scale_args_t sa = {0};
    sa.n = N; sa.mul = 3; sa.add = 5; sa.src = d; sa.dst = out;
    ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_scale, grid1(f.k_scale, N), sa));

    uint32_t host[N];
    ASSERT_TRUE(gpu_read(&f, f.compute, host, out, sizeof(host)));
    for (int i = 0; i < N; ++i) EXPECT_EQ((uint32_t)(i * 3 + 5), host[i]);

    md_gpu_free(f.compute, d);
    md_gpu_free(f.compute, out);
    gpu_close(&f);
}

/* Waiting on a sync from your own stream is a no-op, and a none sync is too. */
UTEST(gpu, wait_on_none_and_self_is_noop) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 64 };
    md_gpu_addr_t d = gpu_alloc(&f, f.compute, N * sizeof(uint32_t));
    ASSERT_TRUE(d != 0);

    md_gpu_stream_wait(f.compute, md_gpu_sync_none());

    fill_args_t a = {0};
    a.n = N; a.base = 42; a.dst = d;
    ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_fill, grid1(f.k_fill, N), a));

    md_gpu_sync_t self = md_gpu_stream_record(f.compute);
    md_gpu_stream_wait(f.compute, self);          /* must not deadlock */

    uint32_t host[N];
    ASSERT_TRUE(gpu_read(&f, f.compute, host, d, sizeof(host)));
    for (int i = 0; i < N; ++i) EXPECT_EQ((uint32_t)(42 + i), host[i]);

    md_gpu_free(f.compute, d);
    gpu_close(&f);
}

/* Concurrency comes from streams. Two streams writing disjoint buffers must
   both land, and each stream is internally ordered. */
UTEST(gpu, two_streams_are_independent) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 512 };
    md_gpu_stream_t sa = md_gpu_stream_create(f.dev, MD_GPU_STREAM_COMPUTE, "a");
    md_gpu_stream_t sb = md_gpu_stream_create(f.dev, MD_GPU_STREAM_COMPUTE, "b");
    ASSERT_TRUE(sa && sb);

    md_gpu_addr_t a = gpu_alloc(&f, sa, N * sizeof(uint32_t));
    md_gpu_addr_t b = gpu_alloc(&f, sb, N * sizeof(uint32_t));
    ASSERT_TRUE(a && b);

    fill_args_t fa = {0}; fa.n = N; fa.base = 10; fa.dst = a;
    fill_args_t fb = {0}; fb.n = N; fb.base = 20; fb.dst = b;
    ASSERT_TRUE(MD_GPU_LAUNCH(sa, f.k_fill, grid1(f.k_fill, N), fa));
    ASSERT_TRUE(MD_GPU_LAUNCH(sb, f.k_fill, grid1(f.k_fill, N), fb));

    scale_args_t da = {0}; da.n = N; da.mul = 2; da.src = a; da.dst = a;
    scale_args_t db = {0}; db.n = N; db.mul = 3; db.src = b; db.dst = b;
    ASSERT_TRUE(MD_GPU_LAUNCH(sa, f.k_scale, grid1(f.k_scale, N), da));
    ASSERT_TRUE(MD_GPU_LAUNCH(sb, f.k_scale, grid1(f.k_scale, N), db));

    uint32_t ha[N], hb[N];
    ASSERT_TRUE(gpu_read(&f, sa, ha, a, sizeof(ha)));
    ASSERT_TRUE(gpu_read(&f, sb, hb, b, sizeof(hb)));
    for (int i = 0; i < N; ++i) {
        EXPECT_EQ((uint32_t)((10 + i) * 2), ha[i]);
        EXPECT_EQ((uint32_t)((20 + i) * 3), hb[i]);
    }

    md_gpu_free(sa, a);
    md_gpu_free(sb, b);
    md_gpu_stream_destroy(sa);
    md_gpu_stream_destroy(sb);
    gpu_close(&f);
}

/* =========================================================================
   Explicit ordering
   ========================================================================= */

/* In EXPLICIT mode md_gpu inserts nothing; the caller's stage barriers carry
   the dependencies. Independent launches need none. */
UTEST(gpu, explicit_ordering_with_stage_barriers) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 1 << 16, ROUNDS = 4 };
    md_gpu_addr_t a = gpu_alloc(&f, f.compute, N * sizeof(uint32_t));
    md_gpu_addr_t b = gpu_alloc(&f, f.compute, N * sizeof(uint32_t));
    md_gpu_addr_t c = gpu_alloc(&f, f.compute, N * sizeof(uint32_t));
    ASSERT_TRUE(a && b && c);
    uint32_t* host = (uint32_t*)malloc(2 * N * sizeof(uint32_t));
    ASSERT_TRUE(host != NULL);

    md_gpu_stream_set_ordering(f.compute, MD_GPU_ORDER_EXPLICIT);
    EXPECT_EQ((int)MD_GPU_ORDER_EXPLICIT, (int)md_gpu_stream_ordering(f.compute));

    for (uint32_t r = 0; r < ROUNDS; ++r) {
        /* Two independent producers: no barrier between them. */
        fill_args_t fa = {0}; fa.n = N; fa.base = 1000 * (r + 1); fa.dst = a;
        fill_args_t fb = {0}; fb.n = N; fb.base = 7 * (r + 1);    fb.dst = b;
        ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_fill, grid1(f.k_fill, N), fa));
        ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_fill, grid1(f.k_fill, N), fb));
        md_gpu_barrier(f.compute, MD_GPU_STAGE_COMPUTE, MD_GPU_STAGE_COMPUTE);

        /* c = a * 2 + 1 consumes a. */
        scale_args_t sa = {0}; sa.n = N; sa.mul = 2; sa.add = 1; sa.src = a; sa.dst = c;
        ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_scale, grid1(f.k_scale, N), sa));
        md_gpu_barrier(f.compute, MD_GPU_STAGE_COMPUTE, MD_GPU_STAGE_TRANSFER);

        md_gpu_mem_t rb = md_gpu_malloc(f.compute, f.read_pool, 2 * N * sizeof(uint32_t));
        ASSERT_TRUE(rb.cpu != NULL);
        ASSERT_TRUE(md_gpu_copy(f.compute, rb.gpu, c, N * sizeof(uint32_t)));
        ASSERT_TRUE(md_gpu_copy(f.compute, rb.gpu + N * sizeof(uint32_t), b, N * sizeof(uint32_t)));
        md_gpu_stream_sync(f.compute);
        memcpy(host, rb.cpu, 2 * N * sizeof(uint32_t));
        md_gpu_free(f.compute, rb.gpu);
        /* The copies above wrote rb; the next round's fills overwrite a, b and
           c, which the copies read -- order them. */
        md_gpu_barrier(f.compute, MD_GPU_STAGE_TRANSFER, MD_GPU_STAGE_COMPUTE);

        for (int i = 0; i < N; i += 131) {
            ASSERT_EQ((1000 * (r + 1) + (uint32_t)i) * 2 + 1, host[i]);
            ASSERT_EQ(7 * (r + 1) + (uint32_t)i, host[N + i]);
        }
    }

    md_gpu_stream_set_ordering(f.compute, MD_GPU_ORDER_IMPLICIT);
    free(host);
    md_gpu_free(f.compute, a);
    md_gpu_free(f.compute, b);
    md_gpu_free(f.compute, c);
    gpu_close(&f);
}

/* Switching back to IMPLICIT orders the next operation after everything in the
   explicit region, with no barrier from the caller. */
UTEST(gpu, explicit_region_end_is_ordered) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 1 << 18, ROUNDS = 4 };
    md_gpu_addr_t a = gpu_alloc(&f, f.compute, N * sizeof(uint32_t));
    md_gpu_addr_t b = gpu_alloc(&f, f.compute, N * sizeof(uint32_t));
    ASSERT_TRUE(a && b);
    uint32_t* host = (uint32_t*)malloc(N * sizeof(uint32_t));
    ASSERT_TRUE(host != NULL);

    for (uint32_t r = 0; r < ROUNDS; ++r) {
        md_gpu_stream_set_ordering(f.compute, MD_GPU_ORDER_EXPLICIT);
        fill_args_t fa = {0}; fa.n = N; fa.base = 100 * (r + 1); fa.dst = a;
        ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_fill, grid1(f.k_fill, N), fa));
        md_gpu_stream_set_ordering(f.compute, MD_GPU_ORDER_IMPLICIT);

        scale_args_t sa = {0}; sa.n = N; sa.mul = 1; sa.add = 5; sa.src = a; sa.dst = b;
        ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_scale, grid1(f.k_scale, N), sa));
        ASSERT_TRUE(gpu_read(&f, f.compute, host, b, N * sizeof(uint32_t)));
        for (int i = 0; i < N; i += 257) ASSERT_EQ(100 * (r + 1) + (uint32_t)i + 5, host[i]);
    }

    free(host);
    md_gpu_free(f.compute, a);
    md_gpu_free(f.compute, b);
    gpu_close(&f);
}

/* =========================================================================
   Textures
   ========================================================================= */

UTEST(gpu, texture_kernel_write_then_read) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { D = 8, VOXELS = D * D * D };
    md_gpu_texture_t tex = gpu_volume(&f, f.compute, D, MD_GPU_TEX_STORAGE);
    ASSERT_TRUE(tex != NULL);

    const md_gpu_texture_desc_t* back = md_gpu_texture_desc(tex);
    ASSERT_TRUE(back != NULL);
    EXPECT_EQ((uint32_t)D, back->width);
    EXPECT_EQ((uint32_t)D, back->depth_or_layers);
    EXPECT_EQ(1u, back->mip_levels);                  /* normalised from 0 */

    tex_args_t ta = {0};
    ta.dim[0] = D; ta.dim[1] = D; ta.dim[2] = D;
    ta.tex = md_gpu_texture_storage(tex, 0);
    ta.scale = 2.0f;
    ASSERT_TRUE(ta.tex.handle != 0);
    ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_tex_write, md_gpu_grid_for(f.k_tex_write, D, D, D), ta));

    float host[VOXELS];
    memset(host, 0, sizeof(host));
    ASSERT_TRUE(gpu_read_tex(&f, f.compute, host, tex, NULL));
    for (int i = 0; i < VOXELS; ++i) EXPECT_EQ((float)i * 2.0f, host[i]);

    md_gpu_texture_destroy(tex);
    md_gpu_device_poll(f.dev);
    gpu_close(&f);
}

UTEST(gpu, texture_upload_then_kernel_read) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { D = 8, VOXELS = D * D * D };
    md_gpu_texture_t tex = gpu_volume(&f, f.compute, D, MD_GPU_TEX_STORAGE);
    ASSERT_TRUE(tex != NULL);

    float src[VOXELS];
    for (int i = 0; i < VOXELS; ++i) src[i] = (float)(VOXELS - i);
    ASSERT_EQ(sizeof(src), md_gpu_texture_region_size(tex, NULL));
    ASSERT_TRUE(md_gpu_upload_texture(f.compute, tex, NULL, src, sizeof(src)));

    md_gpu_addr_t d = gpu_alloc(&f, f.compute, sizeof(src));
    ASSERT_TRUE(d != 0);

    tex_read_args_t ra = {0};
    ra.dim[0] = D; ra.dim[1] = D; ra.dim[2] = D;
    ra.tex = md_gpu_texture_storage(tex, 0);
    ra.dst = d;
    ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_tex_read, md_gpu_grid_for(f.k_tex_read, D, D, D), ra));

    float host[VOXELS];
    memset(host, 0, sizeof(host));
    ASSERT_TRUE(gpu_read(&f, f.compute, host, d, sizeof(host)));
    for (int i = 0; i < VOXELS; ++i) EXPECT_EQ(src[i], host[i]);

    md_gpu_free(f.compute, d);
    md_gpu_texture_destroy(tex);
    md_gpu_device_poll(f.dev);
    gpu_close(&f);
}

/* A texture in a buffer: copy_to_texture / copy_from_texture move data
   between a texture and device memory without touching the host. */
UTEST(gpu, texture_buffer_copies) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { D = 8, VOXELS = D * D * D };
    md_gpu_texture_t tex = gpu_volume(&f, f.compute, D, MD_GPU_TEX_STORAGE);
    ASSERT_TRUE(tex != NULL);
    md_gpu_addr_t a = gpu_alloc(&f, f.compute, VOXELS * sizeof(float));
    md_gpu_addr_t b = gpu_alloc(&f, f.compute, VOXELS * sizeof(float));
    ASSERT_TRUE(a && b);

    float src[VOXELS];
    for (int i = 0; i < VOXELS; ++i) src[i] = 0.25f * (float)i;
    ASSERT_TRUE(md_gpu_upload(f.compute, a, src, sizeof(src)));
    ASSERT_TRUE(md_gpu_copy_to_texture(f.compute, tex, NULL, a));
    ASSERT_TRUE(md_gpu_copy_from_texture(f.compute, b, tex, NULL));

    float host[VOXELS];
    ASSERT_TRUE(gpu_read(&f, f.compute, host, b, sizeof(host)));
    for (int i = 0; i < VOXELS; ++i) EXPECT_EQ(src[i], host[i]);

    md_gpu_free(f.compute, a);
    md_gpu_free(f.compute, b);
    md_gpu_texture_destroy(tex);
    md_gpu_device_poll(f.dev);
    gpu_close(&f);
}

/* Destroying a texture while work that used it is still in flight must be
   safe; its slots come back after a poll once that work has completed. */
UTEST(gpu, texture_destroy_while_in_flight) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { D = 16 };
    for (int iter = 0; iter < 8; ++iter) {
        md_gpu_texture_t tex = gpu_volume(&f, f.compute, D, MD_GPU_TEX_STORAGE);
        ASSERT_TRUE(tex != NULL);
        tex_args_t ta = {0};
        ta.dim[0] = D; ta.dim[1] = D; ta.dim[2] = D;
        ta.tex = md_gpu_texture_storage(tex, 0);
        ta.scale = 1.0f;
        ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_tex_write, md_gpu_grid_for(f.k_tex_write, D, D, D), ta));
        md_gpu_stream_flush(f.compute);
        md_gpu_texture_destroy(tex);       /* still in flight -- legal */
        md_gpu_device_poll(f.dev);
    }

    md_gpu_stream_sync(f.compute);
    md_gpu_device_poll(f.dev);
    gpu_close(&f);
}

/* A 3D texture of depth 1 is still 3D. It used to be silently created as 2D,
   which a RWTexture3D handle then addressed as the wrong view type. */
UTEST(gpu, texture_3d_with_depth_one_stays_3d) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { W = 8, H = 4, TEXELS = W * H };
    md_gpu_texture_desc_t td = {0};
    td.type = MD_GPU_TEX_3D; td.format = MD_GPU_FORMAT_R32_FLOAT; td.usage = MD_GPU_TEX_STORAGE;
    td.width = W; td.height = H; td.depth_or_layers = 1;
    md_gpu_texture_t tex = md_gpu_texture_create(f.compute, f.pool, &td);
    ASSERT_TRUE(tex != NULL);
    EXPECT_EQ((int)MD_GPU_TEX_3D, (int)md_gpu_texture_desc(tex)->type);

    tex_args_t ta = {0};
    ta.dim[0] = W; ta.dim[1] = H; ta.dim[2] = 1;
    ta.tex = md_gpu_texture_storage(tex, 0);
    ta.scale = 3.0f;
    ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_tex_write, md_gpu_grid_for(f.k_tex_write, W, H, 1), ta));

    float host[TEXELS];
    ASSERT_TRUE(gpu_read_tex(&f, f.compute, host, tex, NULL));
    for (int i = 0; i < TEXELS; ++i) EXPECT_EQ((float)i * 3.0f, host[i]);

    md_gpu_texture_destroy(tex);
    gpu_close(&f);
}

/* Invalid descriptions fail at creation, with a reason. */
UTEST(gpu, texture_invalid_descs_are_rejected) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_texture_desc_t td = {0};
    /* Zero-initialised: no type. */
    EXPECT_TRUE(md_gpu_texture_create(f.compute, f.pool, &td) == NULL);
    EXPECT_TRUE(md_gpu_last_error() != NULL);

    /* A 2D texture with depth. */
    td.type = MD_GPU_TEX_2D; td.format = MD_GPU_FORMAT_R32_FLOAT; td.usage = MD_GPU_TEX_STORAGE;
    td.width = 4; td.height = 4; td.depth_or_layers = 4;
    EXPECT_TRUE(md_gpu_texture_create(f.compute, f.pool, &td) == NULL);

    /* No usage. */
    td.depth_or_layers = 1; td.usage = 0;
    EXPECT_TRUE(md_gpu_texture_create(f.compute, f.pool, &td) == NULL);

    /* A host-visible pool. */
    td.usage = MD_GPU_TEX_STORAGE;
    EXPECT_TRUE(md_gpu_texture_create(f.compute, f.read_pool, &td) == NULL);

    /* sRGB cannot be a storage image; the error names the format. */
    td.format = MD_GPU_FORMAT_RGBA8_SRGB;
    EXPECT_TRUE(md_gpu_texture_create(f.compute, f.pool, &td) == NULL);
    const char* err = md_gpu_last_error();
    ASSERT_TRUE(err != NULL);
    EXPECT_TRUE(strstr(err, "RGBA8_SRGB") != NULL);

    gpu_close(&f);
}

/* Each mip level gets its own storage handle; the sampled handle covers all. */
UTEST(gpu, texture_mip_storage_handles) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { D = 8, H = D / 2, MIP1 = H * H * H };
    md_gpu_texture_desc_t td = {0};
    td.type = MD_GPU_TEX_3D; td.format = MD_GPU_FORMAT_R32_FLOAT;
    td.usage = MD_GPU_TEX_STORAGE | MD_GPU_TEX_SAMPLED;
    td.width = D; td.height = D; td.depth_or_layers = D; td.mip_levels = 2;
    md_gpu_texture_t tex = md_gpu_texture_create(f.compute, f.pool, &td);
    ASSERT_TRUE(tex != NULL);

    md_gpu_storage_tex_t m0 = md_gpu_texture_storage(tex, 0);
    md_gpu_storage_tex_t m1 = md_gpu_texture_storage(tex, 1);
    md_gpu_storage_tex_t m2 = md_gpu_texture_storage(tex, 2);
    md_gpu_sampled_tex_t sm = md_gpu_texture_sampled(tex);
    EXPECT_TRUE(m0.handle != 0);
    EXPECT_TRUE(m1.handle != 0);
    EXPECT_TRUE(m0.handle != m1.handle);
    EXPECT_EQ(0u, (unsigned)m2.handle);   /* out of range */
    EXPECT_TRUE(sm.handle != 0);

    /* Write mip 1 through its handle, read it back by region. */
    tex_args_t ta = {0};
    ta.dim[0] = H; ta.dim[1] = H; ta.dim[2] = H;
    ta.tex = m1; ta.scale = 5.0f;
    ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_tex_write, md_gpu_grid_for(f.k_tex_write, H, H, H), ta));

    md_gpu_tex_region_t r = {0};
    r.mip = 1;
    EXPECT_EQ((size_t)(MIP1 * sizeof(float)), md_gpu_texture_region_size(tex, &r));
    float host[MIP1];
    ASSERT_TRUE(gpu_read_tex(&f, f.compute, host, tex, &r));
    for (int i = 0; i < MIP1; ++i) EXPECT_EQ((float)i * 5.0f, host[i]);

    r.mip = 2;
    EXPECT_EQ(0u, (unsigned)md_gpu_texture_region_size(tex, &r));

    md_gpu_texture_destroy(tex);
    gpu_close(&f);
}

/* A sampled handle plus a sampler handle, read through SampleLevel. */
UTEST(gpu, texture_sampled_read_through_sampler) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { D = 8, VOXELS = D * D * D };
    md_gpu_texture_t tex = gpu_volume(&f, f.compute, D, MD_GPU_TEX_SAMPLED);
    ASSERT_TRUE(tex != NULL);
    EXPECT_EQ(0u, (unsigned)md_gpu_texture_storage(tex, 0).handle);   /* not a storage texture */

    float src[VOXELS];
    for (int i = 0; i < VOXELS; ++i) src[i] = (float)(i * 7 % 13);
    ASSERT_TRUE(md_gpu_upload_texture(f.compute, tex, NULL, src, sizeof(src)));

    md_gpu_addr_t d = gpu_alloc(&f, f.compute, sizeof(src));
    ASSERT_TRUE(d != 0);
    sample_args_t sa = {0};
    sa.dim[0] = D; sa.dim[1] = D; sa.dim[2] = D;
    sa.tex = md_gpu_texture_sampled(tex);
    sa.smp = md_gpu_sampler(f.dev, NULL);          /* nearest, clamp */
    sa.dst = d;
    ASSERT_TRUE(sa.tex.handle != 0 && sa.smp.handle != 0);
    ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_sample, md_gpu_grid_for(f.k_sample, D, D, D), sa));

    float host[VOXELS];
    ASSERT_TRUE(gpu_read(&f, f.compute, host, d, sizeof(host)));
    for (int i = 0; i < VOXELS; ++i) EXPECT_EQ(src[i], host[i]);

    md_gpu_free(f.compute, d);
    md_gpu_texture_destroy(tex);
    gpu_close(&f);
}

/* Samplers are cached values: same desc, same handle; nothing to destroy. */
UTEST(gpu, samplers_are_cached) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_sampler_desc_t lin = {0};
    lin.min_filter = lin.mag_filter = lin.mip_filter = MD_GPU_FILTER_LINEAR;
    md_gpu_sampler_t a = md_gpu_sampler(f.dev, &lin);
    md_gpu_sampler_t b = md_gpu_sampler(f.dev, &lin);
    md_gpu_sampler_t c = md_gpu_sampler(f.dev, NULL);
    md_gpu_sampler_desc_t zero = {0};
    md_gpu_sampler_t d = md_gpu_sampler(f.dev, &zero);
    EXPECT_TRUE(a.handle != 0);
    EXPECT_EQ(a.handle, b.handle);
    EXPECT_TRUE(c.handle != 0);
    EXPECT_TRUE(c.handle != a.handle);
    EXPECT_EQ(c.handle, d.handle);                 /* NULL means the zero desc */
    gpu_close(&f);
}

/* Render targets need no heap slot, and depth formats are accepted. */
UTEST(gpu, texture_render_targets) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_texture_desc_t td = {0};
    td.type = MD_GPU_TEX_2D; td.format = MD_GPU_FORMAT_RGBA8_UNORM; td.usage = MD_GPU_TEX_RENDER_TARGET;
    td.width = 64; td.height = 32;
    md_gpu_texture_t color = md_gpu_texture_create(f.compute, f.pool, &td);
    ASSERT_TRUE(color != NULL);
    EXPECT_EQ(0u, (unsigned)md_gpu_texture_storage(color, 0).handle);
    EXPECT_EQ(0u, (unsigned)md_gpu_texture_sampled(color).handle);

    td.format = MD_GPU_FORMAT_D32_FLOAT; td.usage = MD_GPU_TEX_RENDER_TARGET | MD_GPU_TEX_SAMPLED;
    md_gpu_texture_t depth = md_gpu_texture_create(f.compute, f.pool, &td);
    ASSERT_TRUE(depth != NULL);
    EXPECT_TRUE(md_gpu_texture_sampled(depth).handle != 0);

    /* Depth data moves through buffers like any other texture. */
    enum { TEXELS = 64 * 32 };
    float* src = (float*)malloc(TEXELS * sizeof(float));
    float* back = (float*)malloc(TEXELS * sizeof(float));
    ASSERT_TRUE(src && back);
    for (int i = 0; i < TEXELS; ++i) src[i] = (float)i / (float)TEXELS;
    ASSERT_EQ((size_t)(TEXELS * 4), md_gpu_texture_region_size(depth, NULL));
    ASSERT_TRUE(md_gpu_upload_texture(f.compute, depth, NULL, src, TEXELS * sizeof(float)));
    ASSERT_TRUE(gpu_read_tex(&f, f.compute, back, depth, NULL));
    for (int i = 0; i < TEXELS; ++i) EXPECT_EQ(src[i], back[i]);
    free(src); free(back);

    md_gpu_texture_destroy(color);
    md_gpu_texture_destroy(depth);
    gpu_close(&f);
}

/* 2D arrays: the region's z axis addresses layers. */
UTEST(gpu, texture_2d_array_layers) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { W = 8, H = 8, L = 3, LAYER = W * H };
    md_gpu_texture_desc_t td = {0};
    td.type = MD_GPU_TEX_2D_ARRAY; td.format = MD_GPU_FORMAT_R32_UINT; td.usage = MD_GPU_TEX_SAMPLED;
    td.width = W; td.height = H; td.depth_or_layers = L;
    md_gpu_texture_t tex = md_gpu_texture_create(f.compute, f.pool, &td);
    ASSERT_TRUE(tex != NULL);

    uint32_t src[L][LAYER];
    for (int l = 0; l < L; ++l) for (int i = 0; i < LAYER; ++i) src[l][i] = (uint32_t)(l * 1000 + i);
    ASSERT_TRUE(md_gpu_upload_texture(f.compute, tex, NULL, src, sizeof(src)));

    /* Rewrite the middle layer only. */
    uint32_t mid[LAYER];
    for (int i = 0; i < LAYER; ++i) mid[i] = 0xABCD0000u + (uint32_t)i;
    md_gpu_tex_region_t r = {0};
    r.offset[2] = 1; r.extent[2] = 1;
    ASSERT_TRUE(md_gpu_upload_texture(f.compute, tex, &r, mid, sizeof(mid)));

    uint32_t back[L][LAYER];
    ASSERT_TRUE(gpu_read_tex(&f, f.compute, back, tex, NULL));
    for (int i = 0; i < LAYER; ++i) {
        EXPECT_EQ(src[0][i], back[0][i]);
        EXPECT_EQ(mid[i],    back[1][i]);
        EXPECT_EQ(src[2][i], back[2][i]);
    }

    md_gpu_texture_destroy(tex);
    gpu_close(&f);
}

/* Every texture copy so far has covered the whole image, so the region path --
   origin, extent, row and slice pitch -- needs its own test. */
UTEST(gpu, texture_region_subrange_roundtrip) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { D = 8, SUB = 4, VOXELS = D * D * D, SUBVOX = SUB * SUB * SUB };
    md_gpu_texture_t tex = gpu_volume(&f, f.compute, D, MD_GPU_TEX_STORAGE);
    ASSERT_TRUE(tex != NULL);

    float base[VOXELS];
    for (int i = 0; i < VOXELS; ++i) base[i] = (float)i;
    ASSERT_TRUE(md_gpu_upload_texture(f.compute, tex, NULL, base, sizeof(base)));

    /* Overwrite one corner octant. */
    float patch[SUBVOX];
    for (int i = 0; i < SUBVOX; ++i) patch[i] = -1.0f - (float)i;
    md_gpu_tex_region_t region = {0};
    region.offset[0] = SUB; region.offset[1] = SUB; region.offset[2] = SUB;
    region.extent[0] = SUB; region.extent[1] = SUB; region.extent[2] = SUB;
    ASSERT_TRUE(md_gpu_upload_texture(f.compute, tex, &region, patch, sizeof(patch)));

    float back[SUBVOX];
    ASSERT_TRUE(gpu_read_tex(&f, f.compute, back, tex, &region));

    /* The rest of the volume is untouched: sample the opposite corner. */
    md_gpu_tex_region_t other = {0};
    other.extent[0] = SUB; other.extent[1] = SUB; other.extent[2] = SUB;
    float corner[SUBVOX];
    ASSERT_TRUE(gpu_read_tex(&f, f.compute, corner, tex, &other));

    for (int i = 0; i < SUBVOX; ++i) EXPECT_EQ(-1.0f - (float)i, back[i]);
    for (int z = 0; z < SUB; ++z) for (int y = 0; y < SUB; ++y) for (int x = 0; x < SUB; ++x) {
        EXPECT_EQ((float)(z * D * D + y * D + x), corner[z * SUB * SUB + y * SUB + x]);
    }

    md_gpu_texture_destroy(tex);
    gpu_close(&f);
}

/* Copies are sized by the region and checked against the allocation, and a
   region outside the texture is rejected rather than underflowing. */
UTEST(gpu, texture_copies_are_validated) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { D = 8, VOXELS = D * D * D, FULL = VOXELS * sizeof(float) };
    md_gpu_texture_t tex = gpu_volume(&f, f.compute, D, MD_GPU_TEX_STORAGE);
    ASSERT_TRUE(tex != NULL);

    md_gpu_addr_t small = gpu_alloc(&f, f.compute, FULL - sizeof(float));
    md_gpu_addr_t full  = gpu_alloc(&f, f.compute, FULL);
    ASSERT_TRUE(small && full);

    EXPECT_FALSE(md_gpu_copy_from_texture(f.compute, small, tex, NULL));
    EXPECT_FALSE(md_gpu_copy_to_texture(f.compute, tex, NULL, small));
    static float src[VOXELS];
    EXPECT_FALSE(md_gpu_upload_texture(f.compute, tex, NULL, src, FULL - sizeof(float)));

    md_gpu_tex_region_t bad = {0};
    bad.offset[0] = D;                         /* one past the end */
    EXPECT_FALSE(md_gpu_copy_from_texture(f.compute, full, tex, &bad));
    bad.offset[0] = 4; bad.extent[0] = 5;      /* overruns */
    EXPECT_FALSE(md_gpu_copy_from_texture(f.compute, full, tex, &bad));

    /* The exact size works. */
    EXPECT_TRUE(md_gpu_upload_texture(f.compute, tex, NULL, src, FULL));
    EXPECT_TRUE(md_gpu_copy_from_texture(f.compute, full, tex, NULL));
    md_gpu_stream_sync(f.compute);

    md_gpu_free(f.compute, small);
    md_gpu_free(f.compute, full);
    md_gpu_texture_destroy(tex);
    gpu_close(&f);
}

/* A texture handle lives inside the argument struct, and the backends resolve
   it by entirely different means. This separates the two things that can go
   wrong: whether the struct survived the handle (dst[0]) and whether the handle
   itself resolved to a real texture (the readback). */
UTEST(gpu, texture_handle_reaches_shader) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { D = 4, VOXELS = D * D * D, MARKER = 0xC0FFEEu };
    md_gpu_texture_t tex = gpu_volume(&f, f.compute, D, MD_GPU_TEX_STORAGE);
    ASSERT_TRUE(tex != NULL);

    md_gpu_addr_t d = gpu_alloc(&f, f.compute, 2 * sizeof(uint32_t));
    ASSERT_TRUE(d != 0);
    ASSERT_TRUE(md_gpu_memset(f.compute, d, 0, 2 * sizeof(uint32_t)));

    float zero[VOXELS] = {0};
    ASSERT_TRUE(md_gpu_upload_texture(f.compute, tex, NULL, zero, sizeof(zero)));

    tex_probe_args_t a = {0};
    a.n      = 7;
    a.dst    = d;
    a.tex    = md_gpu_texture_storage(tex, 0);
    a.marker = MARKER;
    ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_tex_probe, md_gpu_grid(1, 1, 1), a));

    uint32_t host[2] = {0};
    float voxels[VOXELS];
    ASSERT_TRUE(gpu_read(&f, f.compute, host, d, sizeof(host)));
    ASSERT_TRUE(gpu_read_tex(&f, f.compute, voxels, tex, NULL));

    EXPECT_EQ((uint32_t)MARKER, host[0]);
    EXPECT_EQ(7u, host[1]);
    EXPECT_EQ(1.0f, voxels[0]);
    EXPECT_EQ(2.0f, voxels[1]);

    md_gpu_free(f.compute, d);
    md_gpu_texture_destroy(tex);
    gpu_close(&f);
}

/* A texture destroyed with work outstanding waits on every stream that had
   any -- more streams than any fixed-size list would hold. */
UTEST(gpu, texture_destroy_waits_on_every_stream) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { STREAMS = 12, D = 8 };
    md_gpu_stream_t s[STREAMS];
    for (int i = 0; i < STREAMS; ++i) {
        s[i] = md_gpu_stream_create(f.dev, MD_GPU_STREAM_COMPUTE, "many");
        ASSERT_TRUE(s[i] != NULL);
    }

    md_gpu_texture_t tex = gpu_volume(&f, f.compute, D, MD_GPU_TEX_STORAGE);
    ASSERT_TRUE(tex != NULL);
    md_gpu_sync_t created = md_gpu_stream_record(f.compute);

    for (int i = 0; i < STREAMS; ++i) {
        md_gpu_stream_wait(s[i], created);
        tex_args_t ta = {0};
        ta.dim[0] = D; ta.dim[1] = D; ta.dim[2] = D;
        ta.tex = md_gpu_texture_storage(tex, 0); ta.scale = 1.0f;
        ASSERT_TRUE(MD_GPU_LAUNCH(s[i], f.k_tex_write, md_gpu_grid_for(f.k_tex_write, D, D, D), ta));
        md_gpu_stream_flush(s[i]);
    }

    md_gpu_texture_destroy(tex);
    md_gpu_device_poll(f.dev);

    for (int i = 0; i < STREAMS; ++i) md_gpu_stream_destroy(s[i]);
    md_gpu_device_poll(f.dev);
    gpu_close(&f);
}

/* =========================================================================
   Kernels
   ========================================================================= */

/* The generated descriptor carries [numthreads] and the argument-struct size. */
UTEST(gpu, generated_kernel_descriptors) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_kernel_info_t info = {0};
    ASSERT_TRUE(md_gpu_kernel_info(f.k_fill, &info));
    EXPECT_EQ(64u, info.group_size[0]);
    EXPECT_EQ(1u,  info.group_size[1]);
    EXPECT_EQ(1u,  info.group_size[2]);
    EXPECT_EQ((uint32_t)sizeof(fill_args_t), info.args_size);
    EXPECT_GT(info.max_threads_per_group, 0u);

    ASSERT_TRUE(md_gpu_kernel_info(f.k_tex_write, &info));
    EXPECT_EQ(4u, info.group_size[0]);
    EXPECT_EQ(4u, info.group_size[1]);
    EXPECT_EQ(4u, info.group_size[2]);
    EXPECT_EQ((uint32_t)sizeof(tex_args_t), info.args_size);

    /* Every C mirror in this file agrees with the size the shader build measured. */
    ASSERT_TRUE(md_gpu_kernel_info(f.k_layout, &info));    EXPECT_EQ((uint32_t)sizeof(layout_probe_args_t), info.args_size);
    ASSERT_TRUE(md_gpu_kernel_info(f.k_tex_probe, &info)); EXPECT_EQ((uint32_t)sizeof(tex_probe_args_t), info.args_size);
    ASSERT_TRUE(md_gpu_kernel_info(f.k_sample, &info));    EXPECT_EQ((uint32_t)sizeof(sample_args_t), info.args_size);
    ASSERT_TRUE(md_gpu_kernel_info(f.k_tex_read, &info));  EXPECT_EQ((uint32_t)sizeof(tex_read_args_t), info.args_size);

    md_gpu_grid_t g = md_gpu_grid_for(f.k_tex_write, 9, 8, 1);
    EXPECT_EQ(3u, g.x);
    EXPECT_EQ(2u, g.y);
    EXPECT_EQ(1u, g.z);

    gpu_close(&f);
}

/* A missing or wrong group size is an error on every backend, not a silent
   {1,1,1} that dispatches a fraction of the threads. */
UTEST(gpu, kernel_group_size_is_required_and_checked) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_kernel_desc_t d = md_shader_gpu_test_fill_kernel();
    d.group_size[0] = 0;
    EXPECT_TRUE(md_gpu_kernel_create(f.dev, &d) == NULL);
    EXPECT_TRUE(md_gpu_last_error() != NULL);

#if MD_GPU_BACKEND_VULKAN
    /* Vulkan can also see [numthreads] in the SPIR-V, so a mismatch is caught. */
    d = md_shader_gpu_test_fill_kernel();
    d.group_size[0] = 32;
    EXPECT_TRUE(md_gpu_kernel_create(f.dev, &d) == NULL);
#endif
    gpu_close(&f);
}

/* A launch whose argument block is not the size the kernel was built for is
   rejected before anything is recorded. */
UTEST(gpu, launch_rejects_wrong_argument_size) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_addr_t d = gpu_alloc(&f, f.compute, 64 * sizeof(uint32_t));
    ASSERT_TRUE(d != 0);
    struct { fill_args_t a; uint64_t extra; } big = {{0}, 0};
    big.a.n = 64; big.a.dst = d;
    EXPECT_FALSE(md_gpu_launch(f.compute, f.k_fill, grid1(f.k_fill, 64), &big, sizeof(big)));
    EXPECT_FALSE(md_gpu_launch(f.compute, f.k_fill, grid1(f.k_fill, 64), &big, sizeof(uint32_t)));
    EXPECT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_fill, grid1(f.k_fill, 64), big.a));

    md_gpu_stream_sync(f.compute);
    md_gpu_free(f.compute, d);
    gpu_close(&f);
}

UTEST(gpu, make_grid_and_indirect_launch) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 500 };
    md_gpu_addr_t data  = gpu_alloc(&f, f.compute, N * sizeof(uint32_t));
    md_gpu_addr_t ones  = gpu_alloc(&f, f.compute, N * sizeof(uint32_t));
    md_gpu_addr_t count = gpu_alloc(&f, f.compute, sizeof(uint32_t));
    md_gpu_addr_t grid  = gpu_alloc(&f, f.compute, 3 * sizeof(uint32_t));
    ASSERT_TRUE(data && ones && count && grid);

    fill_args_t fa = {0};
    fa.n = N; fa.base = 0; fa.dst = data;
    ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_fill, grid1(f.k_fill, N), fa));

    scale_args_t za = {0};
    za.n = N; za.mul = 0; za.add = 1; za.src = data; za.dst = ones;
    ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_scale, grid1(f.k_scale, N), za));

    ASSERT_TRUE(md_gpu_memset(f.compute, count, 0, sizeof(uint32_t)));

    sum_args_t sa = {0};
    sa.n = N; sa.src = ones; sa.out_val = count;
    ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_sum, grid1(f.k_sum, N), sa));

    /* The count never touches the CPU: it drives the indirect grid, sized for
       the kernel that will consume it. */
    ASSERT_TRUE(md_gpu_make_grid(f.compute, grid, count, f.k_bump));

    bump_args_t ba = {0};
    ba.n = N; ba.delta = 1000; ba.dst = data;
    ASSERT_TRUE(md_gpu_launch_indirect(f.compute, f.k_bump, grid, &ba, sizeof(ba)));

    uint32_t host_count = 0, host_grid[3] = {0}, host[N];
    ASSERT_TRUE(gpu_read(&f, f.compute, &host_count, count, sizeof(host_count)));
    ASSERT_TRUE(gpu_read(&f, f.compute, host_grid, grid, sizeof(host_grid)));
    ASSERT_TRUE(gpu_read(&f, f.compute, host, data, sizeof(host)));

    EXPECT_EQ((uint32_t)N, host_count);
    EXPECT_EQ((uint32_t)((N + 63) / 64), host_grid[0]);
    EXPECT_EQ(1u, host_grid[1]);
    EXPECT_EQ(1u, host_grid[2]);
    for (int i = 0; i < N; ++i) EXPECT_EQ((uint32_t)(i + 1000), host[i]);

    md_gpu_free(f.compute, data);
    md_gpu_free(f.compute, ones);
    md_gpu_free(f.compute, count);
    md_gpu_free(f.compute, grid);
    gpu_close(&f);
}

/* An indirect grid derived from a device-side count of zero must dispatch
   nothing rather than fault. */
UTEST(gpu, indirect_launch_with_zero_count) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 64 };
    md_gpu_addr_t data  = gpu_alloc(&f, f.compute, N * sizeof(uint32_t));
    md_gpu_addr_t count = gpu_alloc(&f, f.compute, sizeof(uint32_t));
    md_gpu_addr_t grid  = gpu_alloc(&f, f.compute, 3 * sizeof(uint32_t));
    ASSERT_TRUE(data && count && grid);

    fill_args_t fa = {0};
    fa.n = N; fa.base = 0; fa.dst = data;
    ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_fill, grid1(f.k_fill, N), fa));
    ASSERT_TRUE(md_gpu_memset(f.compute, count, 0, sizeof(uint32_t)));
    ASSERT_TRUE(md_gpu_make_grid(f.compute, grid, count, f.k_bump));

    bump_args_t ba = {0};
    ba.n = N; ba.delta = 777; ba.dst = data;
    ASSERT_TRUE(md_gpu_launch_indirect(f.compute, f.k_bump, grid, &ba, sizeof(ba)));

    uint32_t host_grid[3] = {1, 1, 1}, host[N];
    ASSERT_TRUE(gpu_read(&f, f.compute, host_grid, grid, sizeof(host_grid)));
    ASSERT_TRUE(gpu_read(&f, f.compute, host, data, sizeof(host)));

    EXPECT_EQ(0u, host_grid[0]);
    for (int i = 0; i < N; ++i) EXPECT_EQ((uint32_t)i, host[i]);   /* bump did not run */

    md_gpu_free(f.compute, data);
    md_gpu_free(f.compute, count);
    md_gpu_free(f.compute, grid);
    gpu_close(&f);
}

/* =========================================================================
   Host callbacks
   ========================================================================= */

typedef struct {
    int             fired;
    uint32_t        observed;
    const uint32_t* watch;
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
    md_gpu_addr_t d = gpu_alloc(&f, f.compute, N * sizeof(uint32_t));
    md_gpu_mem_t rb = md_gpu_malloc(f.compute, f.read_pool, N * sizeof(uint32_t));
    ASSERT_TRUE(d && rb.cpu);
    memset(rb.cpu, 0, N * sizeof(uint32_t));

    cb_state_t st = {0};
    st.watch = (const uint32_t*)rb.cpu;

    fill_args_t a = {0};
    a.n = N; a.base = 77; a.dst = d;
    ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_fill, grid1(f.k_fill, N), a));
    ASSERT_TRUE(md_gpu_copy(f.compute, rb.gpu, d, N * sizeof(uint32_t)));

    /* Nothing has run the callback yet -- it only fires inside device_poll. */
    ASSERT_TRUE(md_gpu_launch_host_fn(f.compute, gpu_test_callback, &st));
    EXPECT_EQ(0, st.fired);

    md_gpu_stream_sync(f.compute);
    EXPECT_EQ(0, st.fired);       /* still not fired: sync is not poll */

    uint32_t fired = md_gpu_device_poll(f.dev);
    EXPECT_GT(fired, 0u);
    EXPECT_EQ(1, st.fired);
    EXPECT_EQ(77u, st.observed);  /* the copy issued before it is visible */

    md_gpu_device_poll(f.dev);
    EXPECT_EQ(1, st.fired);       /* no re-fire */

    md_gpu_free(f.compute, d);
    md_gpu_free(f.compute, rb.gpu);
    gpu_close(&f);
}

/* A callback observes everything issued into the stream before it, every time.
   Readiness is decided once per poll pass, so a callback can never overtake
   work it was registered after. The loop is what makes a race show. */
UTEST(gpu, host_callback_observes_preceding_copy) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 256, ITERATIONS = 32 };
    md_gpu_addr_t d = gpu_alloc(&f, f.compute, N * sizeof(uint32_t));
    md_gpu_mem_t rb = md_gpu_malloc(f.compute, f.read_pool, N * sizeof(uint32_t));
    ASSERT_TRUE(d && rb.cpu);

    for (uint32_t iter = 0; iter < ITERATIONS; ++iter) {
        fill_args_t a = {0};
        a.n = N; a.base = 1000 + iter; a.dst = d;
        ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_fill, grid1(f.k_fill, N), a));

        cb_state_t st = {0};
        st.watch = (const uint32_t*)rb.cpu;
        ASSERT_TRUE(md_gpu_copy(f.compute, rb.gpu, d, N * sizeof(uint32_t)));
        ASSERT_TRUE(md_gpu_launch_host_fn(f.compute, gpu_test_callback, &st));

        while (!st.fired) md_gpu_device_poll(f.dev);
        EXPECT_EQ(1000u + iter, st.observed);
    }

    md_gpu_free(f.compute, d);
    md_gpu_free(f.compute, rb.gpu);
    gpu_close(&f);
}

UTEST(gpu, sync_on_complete) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 64 };
    md_gpu_addr_t d = gpu_alloc(&f, f.compute, N * sizeof(uint32_t));
    ASSERT_TRUE(d != 0);

    fill_args_t a = {0};
    a.n = N; a.base = 1; a.dst = d;
    ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_fill, grid1(f.k_fill, N), a));

    md_gpu_sync_t sync = md_gpu_stream_record(f.compute);
    cb_state_t st = {0};
    ASSERT_TRUE(md_gpu_sync_on_complete(f.dev, sync, gpu_test_callback, &st));

    /* A none sync fires on the next poll. */
    cb_state_t now = {0};
    ASSERT_TRUE(md_gpu_sync_on_complete(f.dev, md_gpu_sync_none(), gpu_test_callback, &now));

    md_gpu_sync_wait(sync);
    EXPECT_TRUE(md_gpu_sync_is_complete(sync));
    md_gpu_device_poll(f.dev);
    EXPECT_EQ(1, st.fired);
    EXPECT_EQ(1, now.fired);

    md_gpu_free(f.compute, d);
    gpu_close(&f);
}

/* =========================================================================
   Nothing blocks behind a long-running stream

   Async compute jobs here span many frames. Creating or destroying textures,
   kernels and pools on the UI thread must not wait for them -- it used to:
   texture creation submitted a layout transition to compute queue 0 and waited
   on a fence, and pool and kernel destruction idled the whole device.
   ========================================================================= */

/* Launch enough spin work to keep `s` busy for a while. Returns its sync. */
static md_gpu_sync_t gpu_keep_busy(gpu_fixture_t* f, md_gpu_stream_t s, md_gpu_addr_t scratch, uint32_t n, uint32_t iters) {
    spin_args_t sa = {0};
    sa.n = n; sa.iters = iters; sa.dst = scratch;
    for (int i = 0; i < 8; ++i) md_gpu_launch(s, f->k_spin, grid1(f->k_spin, n), &sa, sizeof(sa));
    return md_gpu_stream_record(s);
}

UTEST(gpu, object_lifetime_calls_do_not_wait_for_other_streams) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    /* Lavapipe creates pipelines and image views under the same lock its queue
       thread holds while executing, so on it vkCreateComputePipelines and
       vkCreateImageView themselves wait for running work -- a driver artefact
       that says nothing about md_gpu. Creation is therefore checked on real
       devices only; destruction, pools and allocation everywhere. */
    md_gpu_device_info_t info;
    ASSERT_TRUE(md_gpu_device_info(f.dev, &info));
    const bool software = strstr(info.name, "llvmpipe") != NULL;

    md_gpu_stream_t busy = md_gpu_stream_create(f.dev, MD_GPU_STREAM_COMPUTE, "long job");
    ASSERT_TRUE(busy != NULL);
    enum { N = 4096 };
    md_gpu_addr_t scratch = gpu_alloc(&f, busy, N * sizeof(uint32_t));
    ASSERT_TRUE(scratch != 0);

    /* Objects to destroy while the job runs, created before it starts. */
    md_gpu_pool_t tpool = gpu_pool(f.dev, MD_GPU_MEM_DEVICE, "ui textures");
    ASSERT_TRUE(tpool != NULL);
    md_gpu_texture_desc_t td = {0};
    td.type = MD_GPU_TEX_2D; td.format = MD_GPU_FORMAT_RGBA8_UNORM;
    td.usage = MD_GPU_TEX_SAMPLED | MD_GPU_TEX_STORAGE; td.width = 256; td.height = 256;
    md_gpu_texture_t old_tex = md_gpu_texture_create(f.compute, tpool, &td);
    md_gpu_kernel_desc_t kd = md_shader_gpu_test_fill_kernel();
    md_gpu_kernel_t old_kernel = md_gpu_kernel_create(f.dev, &kd);
    ASSERT_TRUE(old_tex && old_kernel);
    md_gpu_stream_sync(f.compute);

    md_gpu_sync_t job = gpu_keep_busy(&f, busy, scratch, N, 200000);
    if (md_gpu_sync_is_complete(job)) {
        md_gpu_stream_destroy(busy);
        gpu_close(&f);
        UTEST_SKIP("the long-running job finished too quickly for this test to mean anything");
    }

    /* Everything below happens on other streams, or on none. */
    md_gpu_texture_destroy(old_tex);
    md_gpu_kernel_destroy(old_kernel);
    md_gpu_pool_t pool = gpu_pool(f.dev, MD_GPU_MEM_DEVICE, "ui");
    ASSERT_TRUE(pool != NULL);
    ASSERT_TRUE(md_gpu_malloc(f.compute, pool, 1 << 20).gpu != 0);
    md_gpu_pool_destroy(pool);
    md_gpu_pool_destroy(tpool);
    md_gpu_device_poll(f.dev);

    /* If any of those had waited for the device, the job would be done. */
    EXPECT_FALSE(md_gpu_sync_is_complete(job));

    if (!software) {
        md_gpu_pool_t cpool = gpu_pool(f.dev, MD_GPU_MEM_DEVICE, "ui create");
        md_gpu_texture_t tex = md_gpu_texture_create(f.compute, cpool, &td);
        ASSERT_TRUE(tex != NULL);
        md_gpu_kernel_t k = md_gpu_kernel_create(f.dev, &kd);
        ASSERT_TRUE(k != NULL);
        EXPECT_FALSE(md_gpu_sync_is_complete(job));
        md_gpu_kernel_destroy(k);
        md_gpu_texture_destroy(tex);
        md_gpu_pool_destroy(cpool);
    }

    md_gpu_sync_wait(job);
    md_gpu_device_poll(f.dev);
    md_gpu_free(busy, scratch);
    md_gpu_stream_destroy(busy);
    gpu_close(&f);
}

/* Destroying a pool with work in flight is legal and non-blocking; the memory
   goes once that work has completed. */
UTEST(gpu, pool_destroy_with_work_in_flight) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_pool_t pool = gpu_pool(f.dev, MD_GPU_MEM_DEVICE, "in-flight");
    ASSERT_TRUE(pool != NULL);

    enum { N = 65536 };
    md_gpu_addr_t d = md_gpu_malloc(f.compute, pool, N * sizeof(uint32_t)).gpu;
    ASSERT_TRUE(d != 0);
    for (int i = 0; i < 32; ++i) {
        fill_args_t fa = {0};
        fa.n = N; fa.base = (uint32_t)i; fa.dst = d;
        ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_fill, grid1(f.k_fill, N), fa));
    }
    md_gpu_stream_flush(f.compute);

    md_gpu_pool_destroy(pool);        /* must not free memory the GPU is using */

    md_gpu_stream_sync(f.compute);
    md_gpu_device_poll(f.dev);
    gpu_close(&f);
}

/* =========================================================================
   Argument-struct layout

   SPIR-V and MSL lay vectors out differently, and a mismatch does not
   announce itself -- it shows up as a wrong number several launches
   downstream -- so read the fields straight back.
   ========================================================================= */

UTEST(gpu, arg_struct_layout_matches_shader) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 10 };
    md_gpu_addr_t d = gpu_alloc(&f, f.compute, N * sizeof(uint32_t));
    ASSERT_TRUE(d != 0);
    ASSERT_TRUE(md_gpu_memset(f.compute, d, 0xEE, N * sizeof(uint32_t)));

    layout_probe_args_t a = {0};
    a.dim_x = 11; a.dim_y = 22; a.dim_z = 33;
    a.scale = 1.5f;
    a.pair.x = 44; a.pair.y = 55;
    a.v4.x = 2.5f; a.v4.y = 3.5f; a.v4.z = 4.5f; a.v4.w = 5.5f;
    a.dst = d;
    ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_layout, md_gpu_grid(1, 1, 1), a));

    uint32_t host[N];
    ASSERT_TRUE(gpu_read(&f, f.compute, host, d, sizeof(host)));

    float scale, v4[4];
    memcpy(&scale, &host[3], sizeof(scale));
    memcpy(v4, &host[6], sizeof(v4));

    EXPECT_EQ(11u, host[0]);
    EXPECT_EQ(22u, host[1]);
    EXPECT_EQ(33u, host[2]);
    EXPECT_EQ(1.5f, scale);
    EXPECT_EQ(44u, host[4]);
    EXPECT_EQ(55u, host[5]);
    EXPECT_EQ(2.5f, v4[0]);
    EXPECT_EQ(3.5f, v4[1]);
    EXPECT_EQ(4.5f, v4[2]);
    EXPECT_EQ(5.5f, v4[3]);

    md_gpu_free(f.compute, d);
    gpu_close(&f);
}

/* =========================================================================
   Execution hazards

   md_gpu.h's central promise is that every operation observes all writes made
   by the operations before it in the same stream. Neither backend gets that for
   free, and both lose it in a way no small test notices.

     - SIZE. A producer must run long enough for an unordered consumer to get
       ahead of it.
     - ROUNDS. Each round writes a different value, so a stale read reports the
       *previous round's* value rather than zero or garbage.

   Every case is backend-agnostic: nothing here mentions fences, barriers or
   encoders.
   ========================================================================= */

enum {
    HAZ_D      = 128,                        /* 128^3 volume, 8 MiB          */
    HAZ_VOX    = HAZ_D * HAZ_D * HAZ_D,
    HAZ_N      = 1 << 20,                    /* 1 M element buffer, 4 MiB    */
    HAZ_ROUNDS = 4,
    HAZ_STRIDE = 389,                        /* coprime with everything here */
};

/* dispatch -> copy, buffer. */
UTEST(gpu, hazard_dispatch_to_copy_buffer) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_addr_t d = gpu_alloc(&f, f.compute, HAZ_N * sizeof(uint32_t));
    md_gpu_mem_t rb = md_gpu_malloc(f.compute, f.read_pool, HAZ_N * sizeof(uint32_t));
    ASSERT_TRUE(d && rb.cpu);
    const uint32_t* host = (const uint32_t*)rb.cpu;

    for (int r = 0; r < HAZ_ROUNDS; ++r) {
        const uint32_t base = (uint32_t)(r + 1) * 1000000u;
        fill_args_t fa = {0};
        fa.n = HAZ_N; fa.base = base; fa.dst = d;
        ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_fill, grid1(f.k_fill, HAZ_N), fa));
        /* No flush: producer and consumer land in one command buffer. */
        ASSERT_TRUE(md_gpu_copy(f.compute, rb.gpu, d, HAZ_N * sizeof(uint32_t)));
        md_gpu_stream_sync(f.compute);
        for (int i = 0; i < HAZ_N; i += HAZ_STRIDE) ASSERT_EQ(base + (uint32_t)i, host[i]);
    }

    md_gpu_free(f.compute, d);
    md_gpu_free(f.compute, rb.gpu);
    gpu_close(&f);
}

/* dispatch -> copy, texture. The GTO volume readback, reduced. */
UTEST(gpu, hazard_dispatch_to_copy_texture) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_texture_t tex = gpu_volume(&f, f.compute, HAZ_D, MD_GPU_TEX_STORAGE);
    ASSERT_TRUE(tex != NULL);
    md_gpu_mem_t rb = md_gpu_malloc(f.compute, f.read_pool, HAZ_VOX * sizeof(float));
    ASSERT_TRUE(rb.cpu != NULL);
    const float* host = (const float*)rb.cpu;

    for (int r = 0; r < HAZ_ROUNDS; ++r) {
        const float scale = (float)(r + 1);
        tex_args_t ta = {0};
        ta.dim[0] = HAZ_D; ta.dim[1] = HAZ_D; ta.dim[2] = HAZ_D;
        ta.tex = md_gpu_texture_storage(tex, 0); ta.scale = scale;
        ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_tex_write, md_gpu_grid_for(f.k_tex_write, HAZ_D, HAZ_D, HAZ_D), ta));
        ASSERT_TRUE(md_gpu_copy_from_texture(f.compute, rb.gpu, tex, NULL));
        md_gpu_stream_sync(f.compute);
        for (int i = 0; i < HAZ_VOX; i += HAZ_STRIDE) ASSERT_EQ((float)i * scale, host[i]);
    }

    md_gpu_free(f.compute, rb.gpu);
    md_gpu_texture_destroy(tex);
    gpu_close(&f);
}

/* upload -> dispatch, buffer. */
UTEST(gpu, hazard_upload_to_dispatch_buffer) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    uint32_t* src = (uint32_t*)malloc(HAZ_N * sizeof(uint32_t));
    ASSERT_TRUE(src != NULL);
    md_gpu_addr_t a = gpu_alloc(&f, f.compute, HAZ_N * sizeof(uint32_t));
    md_gpu_addr_t b = gpu_alloc(&f, f.compute, HAZ_N * sizeof(uint32_t));
    md_gpu_mem_t rb = md_gpu_malloc(f.compute, f.read_pool, HAZ_N * sizeof(uint32_t));
    ASSERT_TRUE(a && b && rb.cpu);
    const uint32_t* host = (const uint32_t*)rb.cpu;

    for (int r = 0; r < HAZ_ROUNDS; ++r) {
        const uint32_t tag = (uint32_t)(r + 1) * 7u;
        for (int i = 0; i < HAZ_N; ++i) src[i] = (uint32_t)i + tag;
        ASSERT_TRUE(md_gpu_upload(f.compute, a, src, HAZ_N * sizeof(uint32_t)));

        scale_args_t sa = {0};
        sa.n = HAZ_N; sa.mul = 1; sa.add = 0; sa.src = a; sa.dst = b;
        ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_scale, grid1(f.k_scale, HAZ_N), sa));
        ASSERT_TRUE(md_gpu_copy(f.compute, rb.gpu, b, HAZ_N * sizeof(uint32_t)));
        md_gpu_stream_sync(f.compute);
        for (int i = 0; i < HAZ_N; i += HAZ_STRIDE) ASSERT_EQ((uint32_t)i + tag, host[i]);
    }

    md_gpu_free(f.compute, a);
    md_gpu_free(f.compute, b);
    md_gpu_free(f.compute, rb.gpu);
    free(src);
    gpu_close(&f);
}

/* upload -> dispatch, texture. */
UTEST(gpu, hazard_upload_to_dispatch_texture) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { D = 64, VOX = D * D * D };
    md_gpu_texture_t tex = gpu_volume(&f, f.compute, D, MD_GPU_TEX_STORAGE);
    ASSERT_TRUE(tex != NULL);
    float* src = (float*)malloc(VOX * sizeof(float));
    ASSERT_TRUE(src != NULL);
    md_gpu_addr_t d = gpu_alloc(&f, f.compute, VOX * sizeof(float));
    md_gpu_mem_t rb = md_gpu_malloc(f.compute, f.read_pool, VOX * sizeof(float));
    ASSERT_TRUE(d && rb.cpu);
    const float* host = (const float*)rb.cpu;

    for (int r = 0; r < HAZ_ROUNDS; ++r) {
        const float tag = (float)(r + 1) * 1000.0f;
        for (int i = 0; i < VOX; ++i) src[i] = (float)i + tag;
        ASSERT_TRUE(md_gpu_upload_texture(f.compute, tex, NULL, src, VOX * sizeof(float)));

        tex_read_args_t ra = {0};
        ra.dim[0] = D; ra.dim[1] = D; ra.dim[2] = D;
        ra.tex = md_gpu_texture_storage(tex, 0); ra.dst = d;
        ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_tex_read, md_gpu_grid_for(f.k_tex_read, D, D, D), ra));
        ASSERT_TRUE(md_gpu_copy(f.compute, rb.gpu, d, VOX * sizeof(float)));
        md_gpu_stream_sync(f.compute);
        for (int i = 0; i < VOX; i += 61) ASSERT_EQ(src[i], host[i]);
    }

    md_gpu_free(f.compute, d);
    md_gpu_free(f.compute, rb.gpu);
    free(src);
    md_gpu_texture_destroy(tex);
    gpu_close(&f);
}

/* memset -> dispatch. A fill is a transfer operation like any other. */
UTEST(gpu, hazard_memset_to_dispatch) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_addr_t d = gpu_alloc(&f, f.compute, HAZ_N * sizeof(uint32_t));
    md_gpu_mem_t rb = md_gpu_malloc(f.compute, f.read_pool, HAZ_N * sizeof(uint32_t));
    ASSERT_TRUE(d && rb.cpu);
    const uint32_t* host = (const uint32_t*)rb.cpu;

    for (int r = 0; r < HAZ_ROUNDS; ++r) {
        const uint8_t  byte  = (uint8_t)(r + 1);
        const uint32_t word  = (uint32_t)byte * 0x01010101u;
        const uint32_t delta = (uint32_t)(r + 1) * 13u;

        ASSERT_TRUE(md_gpu_memset(f.compute, d, byte, HAZ_N * sizeof(uint32_t)));
        bump_args_t ba = {0};
        ba.n = HAZ_N; ba.delta = delta; ba.dst = d;
        ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_bump, grid1(f.k_bump, HAZ_N), ba));
        ASSERT_TRUE(md_gpu_copy(f.compute, rb.gpu, d, HAZ_N * sizeof(uint32_t)));
        md_gpu_stream_sync(f.compute);
        for (int i = 0; i < HAZ_N; i += HAZ_STRIDE) ASSERT_EQ(word + delta, host[i]);
    }

    md_gpu_free(f.compute, d);
    md_gpu_free(f.compute, rb.gpu);
    gpu_close(&f);
}

/* dispatch -> dispatch through a texture. */
UTEST(gpu, hazard_dispatch_to_dispatch_via_texture) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { D = 64, VOX = D * D * D };
    md_gpu_texture_t tex = gpu_volume(&f, f.compute, D, MD_GPU_TEX_STORAGE);
    ASSERT_TRUE(tex != NULL);
    md_gpu_addr_t d = gpu_alloc(&f, f.compute, VOX * sizeof(float));
    md_gpu_mem_t rb = md_gpu_malloc(f.compute, f.read_pool, VOX * sizeof(float));
    ASSERT_TRUE(d && rb.cpu);
    const float* host = (const float*)rb.cpu;

    for (int r = 0; r < HAZ_ROUNDS; ++r) {
        const float scale = (float)(r + 1) * 0.5f;
        tex_args_t ta = {0};
        ta.dim[0] = D; ta.dim[1] = D; ta.dim[2] = D;
        ta.tex = md_gpu_texture_storage(tex, 0); ta.scale = scale;
        ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_tex_write, md_gpu_grid_for(f.k_tex_write, D, D, D), ta));

        tex_read_args_t ra = {0};
        ra.dim[0] = D; ra.dim[1] = D; ra.dim[2] = D;
        ra.tex = md_gpu_texture_storage(tex, 0); ra.dst = d;
        ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_tex_read, md_gpu_grid_for(f.k_tex_read, D, D, D), ra));

        ASSERT_TRUE(md_gpu_copy(f.compute, rb.gpu, d, VOX * sizeof(float)));
        md_gpu_stream_sync(f.compute);
        for (int i = 0; i < VOX; i += 61) ASSERT_EQ((float)i * scale, host[i]);
    }

    md_gpu_free(f.compute, d);
    md_gpu_free(f.compute, rb.gpu);
    md_gpu_texture_destroy(tex);
    gpu_close(&f);
}

/* Many transitions in one command buffer: transfer, compute, transfer... */
UTEST(gpu, hazard_alternating_chain) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 1 << 18, STEPS = 12 };
    md_gpu_addr_t a = gpu_alloc(&f, f.compute, N * sizeof(uint32_t));
    md_gpu_addr_t b = gpu_alloc(&f, f.compute, N * sizeof(uint32_t));
    ASSERT_TRUE(a && b);

    fill_args_t fa = {0};
    fa.n = N; fa.base = 0; fa.dst = a;
    ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_fill, grid1(f.k_fill, N), fa));
    for (int step = 0; step < STEPS; ++step) {
        ASSERT_TRUE(md_gpu_copy(f.compute, b, a, N * sizeof(uint32_t)));
        scale_args_t sa = {0};
        sa.n = N; sa.mul = 1; sa.add = 1; sa.src = b; sa.dst = a;
        ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_scale, grid1(f.k_scale, N), sa));
    }

    uint32_t* host = (uint32_t*)malloc(N * sizeof(uint32_t));
    ASSERT_TRUE(host != NULL);
    ASSERT_TRUE(gpu_read(&f, f.compute, host, a, N * sizeof(uint32_t)));
    for (int i = 0; i < N; i += 97) ASSERT_EQ((uint32_t)(i + STEPS), host[i]);

    free(host);
    md_gpu_free(f.compute, a);
    md_gpu_free(f.compute, b);
    gpu_close(&f);
}

/* Producer and consumer forced into separate submissions. */
UTEST(gpu, hazard_across_submission_boundary) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_texture_t tex = gpu_volume(&f, f.compute, HAZ_D, MD_GPU_TEX_STORAGE);
    ASSERT_TRUE(tex != NULL);
    md_gpu_mem_t rb = md_gpu_malloc(f.compute, f.read_pool, HAZ_VOX * sizeof(float));
    ASSERT_TRUE(rb.cpu != NULL);
    const float* host = (const float*)rb.cpu;

    for (int r = 0; r < HAZ_ROUNDS; ++r) {
        const float scale = (float)(r + 1) * 3.0f;
        tex_args_t ta = {0};
        ta.dim[0] = HAZ_D; ta.dim[1] = HAZ_D; ta.dim[2] = HAZ_D;
        ta.tex = md_gpu_texture_storage(tex, 0); ta.scale = scale;
        ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_tex_write, md_gpu_grid_for(f.k_tex_write, HAZ_D, HAZ_D, HAZ_D), ta));
        md_gpu_stream_flush(f.compute);     /* producer and consumer now split */
        ASSERT_TRUE(md_gpu_copy_from_texture(f.compute, rb.gpu, tex, NULL));
        md_gpu_stream_sync(f.compute);
        for (int i = 0; i < HAZ_VOX; i += HAZ_STRIDE) ASSERT_EQ((float)i * scale, host[i]);
    }

    md_gpu_free(f.compute, rb.gpu);
    md_gpu_texture_destroy(tex);
    gpu_close(&f);
}

/* Producer and consumer on different streams, joined by a sync. */
UTEST(gpu, hazard_across_streams_via_sync) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_stream_t producer = md_gpu_stream_create(f.dev, MD_GPU_STREAM_COMPUTE, "producer");
    ASSERT_TRUE(producer != NULL);
    md_gpu_texture_t tex = gpu_volume(&f, producer, HAZ_D, MD_GPU_TEX_STORAGE);
    ASSERT_TRUE(tex != NULL);
    md_gpu_mem_t rb = md_gpu_malloc(f.compute, f.read_pool, HAZ_VOX * sizeof(float));
    ASSERT_TRUE(rb.cpu != NULL);
    const float* host = (const float*)rb.cpu;

    for (int r = 0; r < HAZ_ROUNDS; ++r) {
        const float scale = (float)(r + 1) * 7.0f;
        tex_args_t ta = {0};
        ta.dim[0] = HAZ_D; ta.dim[1] = HAZ_D; ta.dim[2] = HAZ_D;
        ta.tex = md_gpu_texture_storage(tex, 0); ta.scale = scale;
        ASSERT_TRUE(MD_GPU_LAUNCH(producer, f.k_tex_write, md_gpu_grid_for(f.k_tex_write, HAZ_D, HAZ_D, HAZ_D), ta));

        md_gpu_sync_t done = md_gpu_stream_record(producer);
        ASSERT_TRUE(md_gpu_sync_is_valid(done));
        md_gpu_stream_wait(f.compute, done);

        ASSERT_TRUE(md_gpu_copy_from_texture(f.compute, rb.gpu, tex, NULL));
        md_gpu_stream_sync(f.compute);
        for (int i = 0; i < HAZ_VOX; i += HAZ_STRIDE) ASSERT_EQ((float)i * scale, host[i]);
        /* The next round's write must not overtake this round's read. */
        md_gpu_stream_wait(producer, md_gpu_stream_record(f.compute));
    }

    md_gpu_free(f.compute, rb.gpu);
    md_gpu_texture_destroy(tex);
    md_gpu_stream_destroy(producer);
    gpu_close(&f);
}

/* =========================================================================
   Synchronisation at scale
   ========================================================================= */

/* Fan more producers into one consumer than any fixed-size wait list or poll
   snapshot would hold. */
UTEST(gpu, wide_cross_stream_fan_in) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { STREAMS = 20, N = 128 };
    md_gpu_stream_t s[STREAMS] = {0};
    md_gpu_addr_t buf[STREAMS] = {0};

    for (int i = 0; i < STREAMS; ++i) {
        s[i] = md_gpu_stream_create(f.dev, MD_GPU_STREAM_COMPUTE, "producer");
        ASSERT_TRUE(s[i] != NULL);
        buf[i] = gpu_alloc(&f, s[i], N * sizeof(uint32_t));
        ASSERT_TRUE(buf[i] != 0);
        fill_args_t fa = {0};
        fa.n = N; fa.base = (uint32_t)(i * 1000); fa.dst = buf[i];
        ASSERT_TRUE(MD_GPU_LAUNCH(s[i], f.k_fill, grid1(f.k_fill, N), fa));
    }

    for (int i = 0; i < STREAMS; ++i) {
        md_gpu_sync_t sync = md_gpu_stream_record(s[i]);
        ASSERT_TRUE(md_gpu_sync_is_valid(sync));
        md_gpu_stream_wait(f.compute, sync);
    }

    md_gpu_mem_t rb = md_gpu_malloc(f.compute, f.read_pool, STREAMS * N * sizeof(uint32_t));
    ASSERT_TRUE(rb.cpu != NULL);
    for (int i = 0; i < STREAMS; ++i) {
        ASSERT_TRUE(md_gpu_copy(f.compute, rb.gpu + (uint64_t)i * N * sizeof(uint32_t), buf[i], N * sizeof(uint32_t)));
    }
    md_gpu_stream_sync(f.compute);

    const uint32_t* host = (const uint32_t*)rb.cpu;
    for (int i = 0; i < STREAMS; ++i) {
        for (int j = 0; j < N; ++j) ASSERT_EQ((uint32_t)(i * 1000 + j), host[i * N + j]);
    }

    md_gpu_free(f.compute, rb.gpu);
    for (int i = 0; i < STREAMS; ++i) {
        md_gpu_free(s[i], buf[i]);
        md_gpu_stream_destroy(s[i]);
    }
    gpu_close(&f);
}

/* Different streams may be driven concurrently from different threads, and
   malloc/free are thread-safe. Every thread owns its stream and its memory, so
   any failure is md_gpu's shared state -- the allocation registry above all. */
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
        md_gpu_addr_t d = md_gpu_malloc(s, c->f->pool, N * sizeof(uint32_t)).gpu;
        if (!d) { snprintf(c->failure, sizeof(c->failure), "allocation: %s", gpu_no_device_reason()); ok = false; break; }

        fill_args_t fa = {0};
        fa.n = N; fa.base = (uint32_t)(c->index * 100000 + it); fa.dst = d;
        if (!MD_GPU_LAUNCH(s, c->f->k_fill, grid1(c->f->k_fill, N), fa)) {
            snprintf(c->failure, sizeof(c->failure), "launch: %s", gpu_no_device_reason()); ok = false; break;
        }
        uint32_t host[N];
        if (!gpu_read(c->f, s, host, d, sizeof(host))) {
            snprintf(c->failure, sizeof(c->failure), "copy: %s", gpu_no_device_reason()); ok = false; break;
        }
        for (int i = 0; i < N && ok; ++i) {
            if (host[i] != (uint32_t)(c->index * 100000 + it + i)) {
                snprintf(c->failure, sizeof(c->failure), "readback data mismatch");
                ok = false;
            }
        }
        md_gpu_free(s, d);
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

/* device_destroy destroys everything created from the device, kernels
   included; the counting allocator says whether that is true. */
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

    md_gpu_pool_t pool = gpu_pool(dev, MD_GPU_MEM_DEVICE, "leaked");
    ASSERT_TRUE(pool != NULL);
    ASSERT_TRUE(md_gpu_malloc(s, pool, 64 * 1024).gpu != 0);

    md_gpu_texture_desc_t td = {0};
    td.type = MD_GPU_TEX_3D; td.format = MD_GPU_FORMAT_R32_FLOAT; td.usage = MD_GPU_TEX_STORAGE;
    td.width = 4; td.height = 4; td.depth_or_layers = 4;
    ASSERT_TRUE(md_gpu_texture_create(s, pool, &td) != NULL);
    ASSERT_TRUE(md_gpu_sampler(dev, NULL).handle != 0);

    md_gpu_kernel_desc_t kd = md_shader_gpu_test_fill_kernel();
    ASSERT_TRUE(md_gpu_kernel_create(dev, &kd) != NULL);

    /* Nothing above is destroyed by hand. */
    md_gpu_device_destroy(dev);
    EXPECT_EQ(0u, (unsigned)stats.live_bytes);
}

/* A pool block records the stream it was freed on so that a later allocation
   can tell whether reuse is safe. Destroying that stream must neither leave
   the record dangling nor strand the block. */
UTEST(gpu, stream_destroy_releases_blocks_it_freed) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_pool_t pool = gpu_pool(f.dev, MD_GPU_MEM_DEVICE, "caching");
    ASSERT_TRUE(pool != NULL);

    enum { N = 1024 };
    md_gpu_stream_t worker = md_gpu_stream_create(f.dev, MD_GPU_STREAM_COMPUTE, "worker");
    ASSERT_TRUE(worker != NULL);

    md_gpu_addr_t a = md_gpu_malloc(worker, pool, N * sizeof(uint32_t)).gpu;
    ASSERT_TRUE(a != 0);
    ASSERT_TRUE(md_gpu_memset(worker, a, 0, N * sizeof(uint32_t)));
    md_gpu_free(worker, a);           /* pending: the memset has not been submitted */

    md_gpu_pool_stats_t before = {0};
    md_gpu_pool_stats(pool, &before);

    md_gpu_stream_destroy(worker);

    md_gpu_addr_t b = md_gpu_malloc(f.compute, pool, N * sizeof(uint32_t)).gpu;
    ASSERT_TRUE(b != 0);
    md_gpu_pool_stats_t after = {0};
    md_gpu_pool_stats(pool, &after);
    EXPECT_EQ(before.bytes_reserved, after.bytes_reserved);   /* the cached block */

    fill_args_t fa = {0};
    fa.n = N; fa.base = 7; fa.dst = b;
    ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_fill, grid1(f.k_fill, N), fa));
    uint32_t host[N];
    ASSERT_TRUE(gpu_read(&f, f.compute, host, b, sizeof(host)));
    for (int i = 0; i < N; ++i) ASSERT_EQ((uint32_t)(7 + i), host[i]);

    md_gpu_free(f.compute, b);
    md_gpu_pool_destroy(pool);
    md_gpu_device_poll(f.dev);
    gpu_close(&f);
}

UTEST(gpu, stress_alloc_launch_free_cycles) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    enum { N = 4096, CYCLES = 64 };
    md_gpu_pool_t pool = gpu_pool(f.dev, MD_GPU_MEM_DEVICE, "stress");
    ASSERT_TRUE(pool != NULL);

    for (int c = 0; c < CYCLES; ++c) {
        md_gpu_addr_t d = md_gpu_malloc(f.compute, pool, N * sizeof(uint32_t)).gpu;
        ASSERT_TRUE(d != 0);
        fill_args_t a = {0};
        a.n = N; a.base = (uint32_t)c; a.dst = d;
        ASSERT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_fill, grid1(f.k_fill, N), a));
        md_gpu_free(f.compute, d);
        if ((c & 7) == 0) md_gpu_device_poll(f.dev);
    }
    md_gpu_stream_sync(f.compute);

    md_gpu_pool_stats_t st = {0};
    md_gpu_pool_stats(pool, &st);
    EXPECT_EQ(0u, (unsigned)st.bytes_in_use);
    EXPECT_LE(st.bytes_reserved, (uint64_t)(N * sizeof(uint32_t) * 4));
    EXPECT_GE(st.reuse_count, st.alloc_count - 1);

    md_gpu_pool_destroy(pool);
    gpu_close(&f);
}

/* =========================================================================
   Edge cases and error reporting
   ========================================================================= */

UTEST(gpu, null_and_zero_arguments_are_tolerated) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    /* Frees and destroys of nothing. */
    md_gpu_free(f.compute, 0);
    md_gpu_texture_destroy(NULL);
    md_gpu_kernel_destroy(NULL);
    md_gpu_pool_destroy(NULL);
    md_gpu_pool_trim(NULL, 0);

    /* Zero-sized transfers are no-ops, not failures. */
    md_gpu_addr_t d = gpu_alloc(&f, f.compute, 256);
    ASSERT_TRUE(d != 0);
    uint32_t host = 0;
    EXPECT_TRUE(md_gpu_upload(f.compute, d, &host, 0));
    EXPECT_TRUE(md_gpu_memset(f.compute, d, 0, 0));
    EXPECT_TRUE(md_gpu_copy(f.compute, d, d, 0));

    /* An empty grid launches nothing and succeeds. */
    fill_args_t fa = {0};
    fa.n = 1; fa.dst = d;
    EXPECT_TRUE(MD_GPU_LAUNCH(f.compute, f.k_fill, md_gpu_grid(0, 1, 1), fa));

    EXPECT_TRUE(md_gpu_texture_desc(NULL) == NULL);
    EXPECT_EQ(0u, (unsigned)md_gpu_texture_storage(NULL, 0).handle);
    EXPECT_EQ(0u, (unsigned)md_gpu_texture_sampled(NULL).handle);

    md_gpu_stream_sync(f.compute);
    md_gpu_free(f.compute, d);
    gpu_close(&f);
}

/* Addresses are checked against live allocations; nothing is inferred. */
UTEST(gpu, addresses_are_validated) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_addr_t d = gpu_alloc(&f, f.compute, 256);
    ASSERT_TRUE(d != 0);
    uint32_t host[128] = {0};

    EXPECT_FALSE(md_gpu_upload(f.compute, d, host, 512));           /* overruns */
    EXPECT_TRUE(md_gpu_last_error() != NULL);
    EXPECT_FALSE(md_gpu_copy(f.compute, d + 128, d, 256));           /* dst overruns */
    EXPECT_FALSE(md_gpu_memset(f.compute, 0x1000, 0, 4));            /* not an allocation */
    EXPECT_FALSE(md_gpu_copy(f.compute, d, (md_gpu_addr_t)(uintptr_t)host, 4));  /* a host pointer is not an address */

    /* A freed address is no longer valid. */
    md_gpu_free(f.compute, d);
    EXPECT_FALSE(md_gpu_memset(f.compute, d, 0, 4));
    /* Freeing requires a stream. */
    md_gpu_addr_t e = gpu_alloc(&f, f.compute, 256);
    ASSERT_TRUE(e != 0);
    md_gpu_free(NULL, e);
    EXPECT_TRUE(md_gpu_last_error() != NULL);
    md_gpu_free(f.compute, e);

    md_gpu_stream_sync(f.compute);
    gpu_close(&f);
}

UTEST(gpu, upload_begin_is_exclusive) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_addr_t d = gpu_alloc(&f, f.compute, 256);
    ASSERT_TRUE(d != 0);
    void* p = md_gpu_upload_begin(f.compute, d, 256);
    ASSERT_TRUE(p != NULL);
    EXPECT_TRUE(md_gpu_upload_begin(f.compute, d, 256) == NULL);   /* one open upload per stream */
    EXPECT_FALSE(md_gpu_memset(f.compute, d, 0, 4));               /* nor other work while it is open */
    EXPECT_TRUE(md_gpu_upload_end(f.compute));
    EXPECT_FALSE(md_gpu_upload_end(f.compute));

    md_gpu_stream_sync(f.compute);
    md_gpu_free(f.compute, d);
    gpu_close(&f);
}

/* A transient page appended to after an earlier, completed submission must not
   be recycled under the data just placed in it. The sizes are chosen around
   the 256 KiB page: the marker lands at ~200 KiB in a page stamped by a
   finished submit; the next upload does not fit; the one after would, from
   offset 0, overrun the marker if that page had been recycled. */
UTEST(gpu, transient_pages_are_not_recycled_under_pending_data) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());
    md_gpu_stream_t s = md_gpu_stream_create(f.dev, MD_GPU_STREAM_COMPUTE, "arena");
    ASSERT_TRUE(s != NULL);

    const size_t big0 = 200 * 1024, big1 = 100 * 1024, big2 = 110 * 1024;
    uint8_t* zeros = (uint8_t*)calloc(1, big2);
    ASSERT_TRUE(zeros != NULL);
    md_gpu_addr_t a = gpu_alloc(&f, s, big0);
    md_gpu_addr_t b = gpu_alloc(&f, s, 64);
    md_gpu_addr_t c = gpu_alloc(&f, s, big1);
    md_gpu_addr_t d = gpu_alloc(&f, s, big2);
    ASSERT_TRUE(a && b && c && d);

    ASSERT_TRUE(md_gpu_upload(s, a, zeros, big0));
    md_gpu_stream_sync(s);                              /* the page is stamped and complete */

    uint8_t marker[64];
    memset(marker, 0xAB, sizeof(marker));
    ASSERT_TRUE(md_gpu_upload(s, b, marker, sizeof(marker)));
    ASSERT_TRUE(md_gpu_upload(s, c, zeros, big1));
    ASSERT_TRUE(md_gpu_upload(s, d, zeros, big2));

    uint8_t got[64];
    ASSERT_TRUE(gpu_read(&f, s, got, b, sizeof(got)));
    EXPECT_EQ(0, memcmp(marker, got, sizeof(got)));

    free(zeros);
    md_gpu_free(s, a); md_gpu_free(s, b); md_gpu_free(s, c); md_gpu_free(s, d);
    md_gpu_stream_destroy(s);
    gpu_close(&f);
}

/* An upload into host-visible memory may take the direct memcpy path only
   when nothing it is ordered after is pending -- including a stream wait. */
UTEST(gpu, direct_upload_respects_stream_wait) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());
    md_gpu_stream_t busy  = md_gpu_stream_create(f.dev, MD_GPU_STREAM_COMPUTE, "busy");
    md_gpu_stream_t other = md_gpu_stream_create(f.dev, MD_GPU_STREAM_COMPUTE, "other");
    ASSERT_TRUE(busy && other);

    enum { N = 256, SPIN_N = 4096 };
    md_gpu_addr_t scratch = gpu_alloc(&f, busy, SPIN_N * sizeof(uint32_t));
    md_gpu_mem_t  x       = md_gpu_malloc(busy, f.write_pool, N * sizeof(uint32_t));
    md_gpu_addr_t y       = gpu_alloc(&f, busy, N * sizeof(uint32_t));
    ASSERT_TRUE(scratch && x.cpu && y);
    md_gpu_stream_sync(busy);
    for (int i = 0; i < N; ++i) ((uint32_t*)x.cpu)[i] = 1;

    /* busy: a long job, then read x into y. other: after that, overwrite x. */
    gpu_keep_busy(&f, busy, scratch, SPIN_N, 20000);
    ASSERT_TRUE(md_gpu_copy(busy, y, x.gpu, N * sizeof(uint32_t)));
    md_gpu_sync_t read_done = md_gpu_stream_record(busy);
    md_gpu_stream_wait(other, read_done);

    uint32_t two[N];
    for (int i = 0; i < N; ++i) two[i] = 2;
    ASSERT_TRUE(md_gpu_upload(other, x.gpu, two, sizeof(two)));
    md_gpu_stream_sync(other);
    md_gpu_stream_sync(busy);

    uint32_t hy[N];
    ASSERT_TRUE(gpu_read(&f, busy, hy, y, sizeof(hy)));
    for (int i = 0; i < N; ++i) EXPECT_EQ(1u, hy[i]);
    for (int i = 0; i < N; ++i) EXPECT_EQ(2u, ((uint32_t*)x.cpu)[i]);

    md_gpu_free(busy, scratch); md_gpu_free(busy, x.gpu); md_gpu_free(busy, y);
    md_gpu_stream_destroy(busy);
    md_gpu_stream_destroy(other);
    gpu_close(&f);
}

UTEST(gpu, device_info_is_sane) {
    gpu_fixture_t f;
    if (!gpu_open(&f)) UTEST_SKIP(gpu_no_device_reason());

    md_gpu_device_info_t info;
    ASSERT_TRUE(md_gpu_device_info(f.dev, &info));
    EXPECT_GT(info.max_threads_per_group, 0u);
    EXPECT_GT(info.preferred_group_multiple, 0u);
    md_gpu_kernel_info_t ki;
    ASSERT_TRUE(md_gpu_kernel_info(f.k_tex_write, &ki));
    EXPECT_LE(ki.group_size[0] * ki.group_size[1] * ki.group_size[2], info.max_threads_per_group);
    EXPECT_EQ(4u, md_gpu_format_texel_size(MD_GPU_FORMAT_R32_FLOAT));
    EXPECT_EQ(16u, md_gpu_format_texel_size(MD_GPU_FORMAT_RGBA32_FLOAT));
    EXPECT_EQ(4u, md_gpu_format_texel_size(MD_GPU_FORMAT_D32_FLOAT_S8_UINT));

    gpu_close(&f);
}

#endif /* MD_ENABLE_GPU */
