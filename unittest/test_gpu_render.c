#include "utest.h"

#include <core/md_gpu.h>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if MD_ENABLE_GPU

#include "gpu_render_test_shaders.inl"
#include "gpu_imgui_test_shaders.inl"

/* Rendering tests: pipelines, passes, dynamic state, draws, presentation.

   Pixel convention used throughout: a W x H target, clip x -> pixel
   (x + 1) / 2 * W, clip y -> pixel (1 - y) / 2 * H. Clip space is +Y up and
   pixel (0, 0) is the top-left texel, on every backend. */

/* Mirrors RArgs in unittest/shaders/gpu_render_test.slang. */
typedef struct {
    md_gpu_float4        color;
    md_gpu_addr_t        pos;
    md_gpu_addr_t        colors;
    md_gpu_addr_t        draws;
    uint32_t             id;
    uint32_t             flags;
    md_gpu_sampled_tex_t tex;
    md_gpu_sampler_t     smp;
    uint32_t             width;
    uint32_t             height;
    uint32_t             _p0;
    uint32_t             _p1;
} rargs_t;

typedef struct { md_gpu_float4 offset; md_gpu_float4 color; } rdraw_t;

enum { RF_VERTEX_COLOR = 1, RF_RECORDS = 2, RF_SAMPLE = 4 };

typedef struct {
    md_gpu_device_t dev;
    md_gpu_stream_t gfx;
    bool            present;
} rfix_t;

static bool r_open(rfix_t* f) {
    memset(f, 0, sizeof(*f));
    md_gpu_device_desc_t dd = {0};
    dd.enable_validation = true;
    dd.label             = "md_gpu render unittest";
    f->dev = md_gpu_device_create(&dd);
    if (!f->dev) return false;
    md_gpu_device_info_t info;
    md_gpu_device_info(f->dev, &info);
    if (!info.supports_graphics) {
        md_gpu_device_destroy(f->dev);
        f->dev = NULL;
        return false;
    }
    f->present = info.supports_present;
    f->gfx = md_gpu_stream_default(f->dev, MD_GPU_STREAM_GRAPHICS);
    return f->gfx != NULL;
}

static const char* r_skip_reason(void) {
    const char* err = md_gpu_last_error();
    return (err && err[0]) ? err : "No GPU device with graphics support";
}

static void r_close(rfix_t* f) {
    md_gpu_device_destroy(f->dev);
}

static md_gpu_texture_t r_target(rfix_t* f, md_gpu_format_t fmt, uint32_t w, uint32_t h, md_gpu_tex_usage_t extra) {
    md_gpu_texture_desc_t td = {0};
    td.type   = MD_GPU_TEX_2D;
    td.format = fmt;
    td.usage  = MD_GPU_TEX_RENDER_TARGET | extra;
    td.width  = w;
    td.height = h;
    td.label  = "render target";
    return md_gpu_texture_create(f->gfx, &td);
}

static md_gpu_pipeline_t r_pipeline(rfix_t* f, md_gpu_shader_t vs, md_gpu_shader_t fs, md_gpu_topology_t topo,
                                    const md_gpu_format_t* colors, uint32_t color_count, md_gpu_format_t depth) {
    md_gpu_pipeline_desc_t pd = {0};
    pd.vertex       = vs;
    pd.fragment     = fs;
    pd.topology     = topo;
    pd.color_count  = color_count;
    for (uint32_t i = 0; i < color_count; ++i) pd.color[i].format = colors[i];
    pd.depth_format = depth;
    pd.label        = "test pipeline";
    return md_gpu_pipeline_create(f->dev, &pd);
}

/* Upload `size` bytes into fresh device memory on the graphics stream. */
static md_gpu_addr_t r_buffer(rfix_t* f, const void* data, size_t size) {
    md_gpu_addr_t a = md_gpu_malloc(f->gfx, MD_GPU_MEM_DEVICE, size).gpu;
    if (a && data) md_gpu_upload(f->gfx, a, data, size);
    return a;
}

/* Whole mip `mip`, layer `layer` of `t`, into `dst`. Synchronises the stream. */
static bool r_read(rfix_t* f, md_gpu_texture_t t, uint32_t mip, uint32_t layer, void* dst, size_t size) {
    md_gpu_tex_region_t r = {0};
    r.mip       = mip;
    r.offset[2] = layer;
    r.extent[2] = 1;
    if (md_gpu_texture_region_size(t, &r) != size) return false;
    md_gpu_mem_t rb = md_gpu_malloc(f->gfx, MD_GPU_MEM_HOST_READ, size);
    if (!rb.cpu) return false;
    bool ok = md_gpu_copy_from_texture(f->gfx, rb.gpu, t, &r);
    md_gpu_stream_sync(f->gfx);
    if (ok) memcpy(dst, rb.cpu, size);
    md_gpu_free(f->gfx, rb.gpu);
    return ok;
}

static uint32_t rgba8(uint8_t r, uint8_t g, uint8_t b, uint8_t a) {
    return (uint32_t)r | ((uint32_t)g << 8) | ((uint32_t)b << 16) | ((uint32_t)a << 24);
}

static bool near8(uint32_t px, uint32_t expect, int tol) {
    for (int c = 0; c < 4; ++c) {
        int a = (int)((px >> (8 * c)) & 0xFF), b = (int)((expect >> (8 * c)) & 0xFF);
        if (abs(a - b) > tol) return false;
    }
    return true;
}

static md_gpu_float4 f4(float x, float y, float z, float w) {
    md_gpu_float4 v; v.x = x; v.y = y; v.z = z; v.w = w; return v;
}

/* A pass with one colour attachment and an optional depth attachment. */
static bool r_begin(rfix_t* f, md_gpu_texture_t color, md_gpu_load_t load, float r, float g, float b, float a,
                    md_gpu_texture_t depth, float clear_depth) {
    md_gpu_render_desc_t rd = {0};
    rd.color_count = color ? 1 : 0;
    rd.color[0].texture = color;
    rd.color[0].load    = load;
    rd.color[0].clear.f32[0] = r; rd.color[0].clear.f32[1] = g;
    rd.color[0].clear.f32[2] = b; rd.color[0].clear.f32[3] = a;
    rd.depth.texture     = depth;
    rd.depth.load        = MD_GPU_LOAD_CLEAR;
    rd.depth.clear_depth = clear_depth;
    rd.label = "test pass";
    return md_gpu_render_begin(f->gfx, &rd);
}

#define W 16
#define H 16

static const md_gpu_float4 full_tri[3] = {
    {-1, -1, 0.5f, 1}, {3, -1, 0.5f, 1}, {-1, 3, 0.5f, 1},
};

/* =========================================================================
   Passes, clears, orientation
   ========================================================================= */

UTEST(gpu_render, clear_colour_and_integer_targets) {
    rfix_t f;
    if (!r_open(&f)) UTEST_SKIP(r_skip_reason());

    md_gpu_texture_t c = r_target(&f, MD_GPU_FORMAT_RGBA8_UNORM, W, H, 0);
    md_gpu_texture_t u = r_target(&f, MD_GPU_FORMAT_R32_UINT, W, H, 0);
    ASSERT_TRUE(c && u);

    md_gpu_render_desc_t rd = {0};
    rd.color_count = 2;
    rd.color[0].texture = c;
    rd.color[0].load    = MD_GPU_LOAD_CLEAR;
    rd.color[0].clear.f32[0] = 1.0f; rd.color[0].clear.f32[3] = 1.0f;
    rd.color[1].texture = u;
    rd.color[1].load    = MD_GPU_LOAD_CLEAR;
    rd.color[1].clear.u32[0] = 0xDEADBEEFu;
    ASSERT_TRUE(md_gpu_render_begin(f.gfx, &rd));
    ASSERT_TRUE(md_gpu_render_end(f.gfx));

    static uint32_t px[W * H], ids[W * H];
    ASSERT_TRUE(r_read(&f, c, 0, 0, px, sizeof(px)));
    ASSERT_TRUE(r_read(&f, u, 0, 0, ids, sizeof(ids)));
    for (int i = 0; i < W * H; ++i) {
        EXPECT_EQ(rgba8(255, 0, 0, 255), px[i]);
        EXPECT_EQ(0xDEADBEEFu, ids[i]);
    }
    md_gpu_texture_destroy(c);
    md_gpu_texture_destroy(u);
    r_close(&f);
}

UTEST(gpu_render, clip_space_is_y_up_and_ccw_is_front) {
    rfix_t f;
    if (!r_open(&f)) UTEST_SKIP(r_skip_reason());

    md_gpu_texture_t c = r_target(&f, MD_GPU_FORMAT_RGBA8_UNORM, W, H, 0);
    const md_gpu_format_t fmt = MD_GPU_FORMAT_RGBA8_UNORM;
    md_gpu_pipeline_t p = r_pipeline(&f, md_shader_gpu_render_test_basic_vs_shader(), md_shader_gpu_render_test_color_fs_shader(),
                                     MD_GPU_TOPOLOGY_TRIANGLES, &fmt, 1, MD_GPU_FORMAT_INVALID);
    ASSERT_TRUE(c && p);

    /* Top-left triangle, counter-clockwise in clip space; bottom-right
       triangle, clockwise. */
    const md_gpu_float4 pos[6] = {
        {-1, 1, 0, 1}, {-1, 0, 0, 1}, {0, 1, 0, 1},
        { 1,-1, 0, 1}, { 0,-1, 0, 1}, {1, 0, 0, 1},
    };
    md_gpu_addr_t vb = r_buffer(&f, pos, sizeof(pos));
    ASSERT_TRUE(vb != 0);

    rargs_t a = {0};
    a.color = f4(0, 1, 0, 1);
    a.pos   = vb;

    struct { md_gpu_cull_t cull; bool cw; bool tl; bool br; } cases[] = {
        {MD_GPU_CULL_NONE,  false, true,  true },
        {MD_GPU_CULL_BACK,  false, true,  false},
        {MD_GPU_CULL_FRONT, false, false, true },
        {MD_GPU_CULL_BACK,  true,  false, true },
    };
    static uint32_t px[W * H];
    for (size_t k = 0; k < sizeof(cases) / sizeof(cases[0]); ++k) {
        ASSERT_TRUE(r_begin(&f, c, MD_GPU_LOAD_CLEAR, 0, 0, 0, 1, NULL, 0));
        md_gpu_draw_state_t ds = {0};
        ds.cull            = cases[k].cull;
        ds.front_clockwise = cases[k].cw;
        md_gpu_set_draw_state(f.gfx, &ds);
        EXPECT_TRUE(MD_GPU_DRAW(f.gfx, p, 6, 1, a));
        ASSERT_TRUE(md_gpu_render_end(f.gfx));
        ASSERT_TRUE(r_read(&f, c, 0, 0, px, sizeof(px)));
        const uint32_t green = rgba8(0, 255, 0, 255), black = rgba8(0, 0, 0, 255);
        EXPECT_EQ(cases[k].tl ? green : black, px[1 * W + 1]);       /* top-left     */
        EXPECT_EQ(cases[k].br ? green : black, px[14 * W + 14]);     /* bottom-right */
        EXPECT_EQ(black, px[1 * W + 14]);                            /* top-right    */
        EXPECT_EQ(black, px[14 * W + 1]);                            /* bottom-left  */
    }
    md_gpu_free(f.gfx, vb);
    md_gpu_pipeline_destroy(p);
    md_gpu_texture_destroy(c);
    r_close(&f);
}

/* =========================================================================
   Depth
   ========================================================================= */

UTEST(gpu_render, depth_test_and_write) {
    rfix_t f;
    if (!r_open(&f)) UTEST_SKIP(r_skip_reason());

    md_gpu_texture_t c = r_target(&f, MD_GPU_FORMAT_RGBA8_UNORM, W, H, 0);
    md_gpu_texture_t d = r_target(&f, MD_GPU_FORMAT_D32_FLOAT, W, H, 0);
    const md_gpu_format_t fmt = MD_GPU_FORMAT_RGBA8_UNORM;
    md_gpu_pipeline_t p = r_pipeline(&f, md_shader_gpu_render_test_basic_vs_shader(), md_shader_gpu_render_test_color_fs_shader(),
                                     MD_GPU_TOPOLOGY_TRIANGLES, &fmt, 1, MD_GPU_FORMAT_D32_FLOAT);
    ASSERT_TRUE(c && d && p);

    /* Four full-screen triangles at depths 0.5, 0.7, 0.3, 0.1. */
    md_gpu_float4 pos[12];
    const float z[4] = {0.5f, 0.7f, 0.3f, 0.1f};
    for (int t = 0; t < 4; ++t) for (int v = 0; v < 3; ++v) { pos[t * 3 + v] = full_tri[v]; pos[t * 3 + v].z = z[t]; }
    md_gpu_addr_t vb = r_buffer(&f, pos, sizeof(pos));

    ASSERT_TRUE(r_begin(&f, c, MD_GPU_LOAD_CLEAR, 0, 0, 0, 1, d, 1.0f));
    md_gpu_draw_state_t ds = {0};
    ds.depth_compare = MD_GPU_COMPARE_LESS;
    ds.depth_write   = true;
    md_gpu_set_draw_state(f.gfx, &ds);
    rargs_t a = {0};
    a.pos = vb;
    a.color = f4(1, 0, 0, 1); EXPECT_TRUE(MD_GPU_DRAW(f.gfx, p, 3, 1, a));      /* 0.5: passes     */
    a.pos = vb + 3 * sizeof(md_gpu_float4);
    a.color = f4(0, 1, 0, 1); EXPECT_TRUE(MD_GPU_DRAW(f.gfx, p, 3, 1, a));      /* 0.7: fails      */
    a.pos = vb + 6 * sizeof(md_gpu_float4);
    a.color = f4(0, 0, 1, 1); EXPECT_TRUE(MD_GPU_DRAW(f.gfx, p, 3, 1, a));      /* 0.3: passes     */
    ds.depth_write = false;
    md_gpu_set_draw_state(f.gfx, &ds);
    a.pos = vb + 9 * sizeof(md_gpu_float4);
    a.color = f4(1, 1, 1, 1); EXPECT_TRUE(MD_GPU_DRAW(f.gfx, p, 3, 1, a));      /* 0.1: colour only */
    ASSERT_TRUE(md_gpu_render_end(f.gfx));

    static uint32_t px[W * H];
    static float    dp[W * H];
    ASSERT_TRUE(r_read(&f, c, 0, 0, px, sizeof(px)));
    ASSERT_TRUE(r_read(&f, d, 0, 0, dp, sizeof(dp)));
    EXPECT_EQ(rgba8(255, 255, 255, 255), px[5 * W + 7]);
    EXPECT_NEAR(0.3f, dp[5 * W + 7], 1e-6f);
    EXPECT_NEAR(0.3f, dp[15 * W + 15], 1e-6f);

    md_gpu_free(f.gfx, vb);
    md_gpu_pipeline_destroy(p);
    md_gpu_texture_destroy(c);
    md_gpu_texture_destroy(d);
    r_close(&f);
}

UTEST(gpu_render, fragment_shader_writes_depth) {
    rfix_t f;
    if (!r_open(&f)) UTEST_SKIP(r_skip_reason());

    md_gpu_texture_t c = r_target(&f, MD_GPU_FORMAT_RGBA8_UNORM, W, H, 0);
    md_gpu_texture_t d = r_target(&f, MD_GPU_FORMAT_D32_FLOAT, W, H, MD_GPU_TEX_SAMPLED);
    const md_gpu_format_t fmt = MD_GPU_FORMAT_RGBA8_UNORM;
    md_gpu_pipeline_t p = r_pipeline(&f, md_shader_gpu_render_test_basic_vs_shader(), md_shader_gpu_render_test_depth_fs_shader(),
                                     MD_GPU_TOPOLOGY_TRIANGLES, &fmt, 1, MD_GPU_FORMAT_D32_FLOAT);
    ASSERT_TRUE(c && d && p);
    md_gpu_addr_t vb = r_buffer(&f, full_tri, sizeof(full_tri));

    ASSERT_TRUE(r_begin(&f, c, MD_GPU_LOAD_CLEAR, 0, 0, 0, 1, d, 1.0f));
    md_gpu_draw_state_t ds = {0};
    ds.depth_compare = MD_GPU_COMPARE_LESS_EQUAL;
    ds.depth_write   = true;
    md_gpu_set_draw_state(f.gfx, &ds);
    rargs_t a = {0};
    a.pos   = vb;
    a.color = f4(0.25f, 0, 0, 1);       /* .x is the depth written */
    EXPECT_TRUE(MD_GPU_DRAW(f.gfx, p, 3, 1, a));
    ASSERT_TRUE(md_gpu_render_end(f.gfx));

    static float dp[W * H];
    ASSERT_TRUE(r_read(&f, d, 0, 0, dp, sizeof(dp)));
    for (int i = 0; i < W * H; ++i) EXPECT_NEAR(0.25f, dp[i], 1e-6f);

    md_gpu_free(f.gfx, vb);
    md_gpu_pipeline_destroy(p);
    md_gpu_texture_destroy(c);
    md_gpu_texture_destroy(d);
    r_close(&f);
}

UTEST(gpu_render, depth_only_pass_without_fragment_shader) {
    rfix_t f;
    if (!r_open(&f)) UTEST_SKIP(r_skip_reason());

    md_gpu_texture_t d = r_target(&f, MD_GPU_FORMAT_D32_FLOAT, W, H, 0);
    md_gpu_shader_t none = {0};
    md_gpu_pipeline_t p = r_pipeline(&f, md_shader_gpu_render_test_basic_vs_shader(), none,
                                     MD_GPU_TOPOLOGY_TRIANGLES, NULL, 0, MD_GPU_FORMAT_D32_FLOAT);
    ASSERT_TRUE(d && p);
    md_gpu_float4 pos[3];
    for (int v = 0; v < 3; ++v) { pos[v] = full_tri[v]; pos[v].z = 0.625f; }
    md_gpu_addr_t vb = r_buffer(&f, pos, sizeof(pos));

    ASSERT_TRUE(r_begin(&f, NULL, MD_GPU_LOAD_CLEAR, 0, 0, 0, 0, d, 1.0f));
    md_gpu_draw_state_t ds = {0};
    ds.depth_compare = MD_GPU_COMPARE_LESS;
    ds.depth_write   = true;
    md_gpu_set_draw_state(f.gfx, &ds);
    rargs_t a = {0};
    a.pos = vb;
    EXPECT_TRUE(MD_GPU_DRAW(f.gfx, p, 3, 1, a));
    ASSERT_TRUE(md_gpu_render_end(f.gfx));

    static float dp[W * H];
    ASSERT_TRUE(r_read(&f, d, 0, 0, dp, sizeof(dp)));
    EXPECT_NEAR(0.625f, dp[0], 1e-6f);
    EXPECT_NEAR(0.625f, dp[W * H - 1], 1e-6f);

    md_gpu_free(f.gfx, vb);
    md_gpu_pipeline_destroy(p);
    md_gpu_texture_destroy(d);
    r_close(&f);
}

UTEST(gpu_render, depth_stencil_format_as_depth_attachment) {
    rfix_t f;
    if (!r_open(&f)) UTEST_SKIP(r_skip_reason());
    md_gpu_texture_t d = r_target(&f, MD_GPU_FORMAT_D32_FLOAT_S8_UINT, W, H, 0);
    md_gpu_shader_t none = {0};
    md_gpu_pipeline_t p = r_pipeline(&f, md_shader_gpu_render_test_basic_vs_shader(), none,
                                     MD_GPU_TOPOLOGY_TRIANGLES, NULL, 0, MD_GPU_FORMAT_D32_FLOAT_S8_UINT);
    ASSERT_TRUE(d && p);
    md_gpu_float4 pos[3];
    for (int v = 0; v < 3; ++v) { pos[v] = full_tri[v]; pos[v].z = 0.375f; }
    md_gpu_addr_t vb = r_buffer(&f, pos, sizeof(pos));
    ASSERT_TRUE(r_begin(&f, NULL, MD_GPU_LOAD_CLEAR, 0, 0, 0, 0, d, 1.0f));
    md_gpu_draw_state_t ds = {0};
    ds.depth_compare = MD_GPU_COMPARE_LESS;
    ds.depth_write   = true;
    md_gpu_set_draw_state(f.gfx, &ds);
    rargs_t a = {0};
    a.pos = vb;
    EXPECT_TRUE(MD_GPU_DRAW(f.gfx, p, 3, 1, a));
    ASSERT_TRUE(md_gpu_render_end(f.gfx));
    static float dp[W * H];
    ASSERT_TRUE(r_read(&f, d, 0, 0, dp, sizeof(dp)));
    EXPECT_NEAR(0.375f, dp[W * H / 2], 1e-6f);
    md_gpu_free(f.gfx, vb);
    md_gpu_pipeline_destroy(p);
    md_gpu_texture_destroy(d);
    r_close(&f);
}

/* A stream destroyed with a pass still open closes it; so does the device. */
UTEST(gpu_render, teardown_with_an_open_pass) {
    rfix_t f;
    if (!r_open(&f)) UTEST_SKIP(r_skip_reason());
    md_gpu_stream_t g2 = md_gpu_stream_create(f.dev, MD_GPU_STREAM_GRAPHICS, "second graphics");
    md_gpu_texture_t c = r_target(&f, MD_GPU_FORMAT_RGBA8_UNORM, W, H, 0);
    ASSERT_TRUE(g2 && c);
    md_gpu_stream_wait(g2, md_gpu_stream_record(f.gfx));
    md_gpu_render_desc_t rd = {0};
    rd.color_count = 1;
    rd.color[0].texture = c;
    rd.color[0].load    = MD_GPU_LOAD_CLEAR;
    ASSERT_TRUE(md_gpu_render_begin(g2, &rd));
    md_gpu_stream_destroy(g2);
    ASSERT_TRUE(md_gpu_render_begin(f.gfx, &rd));
    r_close(&f);
}

/* =========================================================================
   Multiple targets, instancing, blending, viewport and scissor
   ========================================================================= */

UTEST(gpu_render, mrt_with_integer_target_and_instancing) {
    rfix_t f;
    if (!r_open(&f)) UTEST_SKIP(r_skip_reason());

    md_gpu_texture_t c = r_target(&f, MD_GPU_FORMAT_RGBA8_UNORM, W, H, 0);
    md_gpu_texture_t u = r_target(&f, MD_GPU_FORMAT_R32_UINT, W, H, 0);
    const md_gpu_format_t fmts[2] = {MD_GPU_FORMAT_RGBA8_UNORM, MD_GPU_FORMAT_R32_UINT};
    md_gpu_pipeline_t p = r_pipeline(&f, md_shader_gpu_render_test_basic_vs_shader(), md_shader_gpu_render_test_mrt_fs_shader(),
                                     MD_GPU_TOPOLOGY_TRIANGLES, fmts, 2, MD_GPU_FORMAT_INVALID);
    ASSERT_TRUE(c && u && p);

    const md_gpu_float4 pos[3] = {{-1, 1, 0, 1}, {-1, 0, 0, 1}, {0, 1, 0, 1}};
    rdraw_t rec[3];
    for (int i = 0; i < 3; ++i) {
        rec[i].offset = f4(0.5f * (float)i, 0, 0, 0);
        rec[i].color  = f4(i == 0 ? 1.0f : 0.0f, i == 1 ? 1.0f : 0.0f, i == 2 ? 1.0f : 0.0f, 1);
    }
    md_gpu_addr_t vb = r_buffer(&f, pos, sizeof(pos));
    md_gpu_addr_t rb = r_buffer(&f, rec, sizeof(rec));

    md_gpu_render_desc_t rd = {0};
    rd.color_count = 2;
    rd.color[0].texture = c; rd.color[0].load = MD_GPU_LOAD_CLEAR;
    rd.color[1].texture = u; rd.color[1].load = MD_GPU_LOAD_CLEAR; rd.color[1].clear.u32[0] = 0xFFFFFFFFu;
    ASSERT_TRUE(md_gpu_render_begin(f.gfx, &rd));
    rargs_t a = {0};
    a.color = f4(1, 1, 1, 1);
    a.pos   = vb;
    a.draws = rb;
    a.flags = RF_RECORDS;
    a.id    = 100;
    EXPECT_TRUE(MD_GPU_DRAW(f.gfx, p, 3, 3, a));
    ASSERT_TRUE(md_gpu_render_end(f.gfx));

    static uint32_t px[W * H], ids[W * H];
    ASSERT_TRUE(r_read(&f, c, 0, 0, px, sizeof(px)));
    ASSERT_TRUE(r_read(&f, u, 0, 0, ids, sizeof(ids)));
    /* Later instances cover earlier ones where they overlap. */
    EXPECT_EQ(100u, ids[1 * W + 1]);
    EXPECT_EQ(101u, ids[1 * W + 5]);
    EXPECT_EQ(102u, ids[1 * W + 9]);
    EXPECT_EQ(0xFFFFFFFFu, ids[14 * W + 1]);
    EXPECT_EQ(rgba8(255, 0, 0, 255), px[1 * W + 1]);
    EXPECT_EQ(rgba8(0, 255, 0, 255), px[1 * W + 5]);
    EXPECT_EQ(rgba8(0, 0, 255, 255), px[1 * W + 9]);

    md_gpu_free(f.gfx, vb);
    md_gpu_free(f.gfx, rb);
    md_gpu_pipeline_destroy(p);
    md_gpu_texture_destroy(c);
    md_gpu_texture_destroy(u);
    r_close(&f);
}

UTEST(gpu_render, blending_and_blend_constant) {
    rfix_t f;
    if (!r_open(&f)) UTEST_SKIP(r_skip_reason());

    md_gpu_texture_t c = r_target(&f, MD_GPU_FORMAT_RGBA8_UNORM, W, H, 0);
    md_gpu_pipeline_desc_t pd = {0};
    pd.vertex   = md_shader_gpu_render_test_basic_vs_shader();
    pd.fragment = md_shader_gpu_render_test_color_fs_shader();
    pd.color_count = 1;
    pd.color[0].format = MD_GPU_FORMAT_RGBA8_UNORM;
    pd.color[0].blend  = (md_gpu_blend_t){ .enable = true,
        .src_color = MD_GPU_BLEND_SRC_ALPHA, .dst_color = MD_GPU_BLEND_ONE_MINUS_SRC_ALPHA,
        .src_alpha = MD_GPU_BLEND_ONE,       .dst_alpha = MD_GPU_BLEND_ONE_MINUS_SRC_ALPHA };
    md_gpu_pipeline_t alpha = md_gpu_pipeline_create(f.dev, &pd);
    pd.color[0].blend  = (md_gpu_blend_t){ .enable = true,
        .src_color = MD_GPU_BLEND_CONSTANT, .dst_color = MD_GPU_BLEND_ZERO,
        .src_alpha = MD_GPU_BLEND_ONE,      .dst_alpha = MD_GPU_BLEND_ZERO };
    pd.color[0].write_disable = MD_GPU_COLOR_A;
    md_gpu_pipeline_t konst = md_gpu_pipeline_create(f.dev, &pd);
    ASSERT_TRUE(c && alpha && konst);
    md_gpu_addr_t vb = r_buffer(&f, full_tri, sizeof(full_tri));

    static uint32_t px[W * H];
    rargs_t a = {0};
    a.pos = vb;

    ASSERT_TRUE(r_begin(&f, c, MD_GPU_LOAD_CLEAR, 0, 0, 0, 1, NULL, 0));
    a.color = f4(1, 1, 1, 0.5f);
    EXPECT_TRUE(MD_GPU_DRAW(f.gfx, alpha, 3, 1, a));
    ASSERT_TRUE(md_gpu_render_end(f.gfx));
    ASSERT_TRUE(r_read(&f, c, 0, 0, px, sizeof(px)));
    EXPECT_TRUE(near8(px[7 * W + 7], rgba8(128, 128, 128, 255), 1));

    ASSERT_TRUE(r_begin(&f, c, MD_GPU_LOAD_CLEAR, 0, 0, 0, 0.5f, NULL, 0));
    md_gpu_draw_state_t ds = {0};
    ds.blend_constant[0] = 0.25f; ds.blend_constant[1] = 0.5f; ds.blend_constant[2] = 0.75f; ds.blend_constant[3] = 1.0f;
    md_gpu_set_draw_state(f.gfx, &ds);
    a.color = f4(1, 1, 1, 1);
    EXPECT_TRUE(MD_GPU_DRAW(f.gfx, konst, 3, 1, a));
    ASSERT_TRUE(md_gpu_render_end(f.gfx));
    ASSERT_TRUE(r_read(&f, c, 0, 0, px, sizeof(px)));
    /* Colour = constant * source; alpha channel masked, keeps the clear. */
    EXPECT_TRUE(near8(px[7 * W + 7], rgba8(64, 128, 191, 128), 1));

    md_gpu_free(f.gfx, vb);
    md_gpu_pipeline_destroy(alpha);
    md_gpu_pipeline_destroy(konst);
    md_gpu_texture_destroy(c);
    r_close(&f);
}

UTEST(gpu_render, viewport_and_scissor) {
    rfix_t f;
    if (!r_open(&f)) UTEST_SKIP(r_skip_reason());

    md_gpu_texture_t c = r_target(&f, MD_GPU_FORMAT_RGBA8_UNORM, W, H, 0);
    const md_gpu_format_t fmt = MD_GPU_FORMAT_RGBA8_UNORM;
    md_gpu_pipeline_t p = r_pipeline(&f, md_shader_gpu_render_test_basic_vs_shader(), md_shader_gpu_render_test_color_fs_shader(),
                                     MD_GPU_TOPOLOGY_TRIANGLES, &fmt, 1, MD_GPU_FORMAT_INVALID);
    /* A full-screen triangle, then a top-left one. */
    md_gpu_float4 pos[6] = {
        {-1, -1, 0, 1}, {3, -1, 0, 1}, {-1, 3, 0, 1},
        {-1,  1, 0, 1}, {-1, 0, 0, 1}, { 0, 1, 0, 1},
    };
    md_gpu_addr_t vb = r_buffer(&f, pos, sizeof(pos));
    ASSERT_TRUE(c && p && vb);
    static uint32_t px[W * H];
    const uint32_t red = rgba8(255, 0, 0, 255), black = rgba8(0, 0, 0, 255);

    /* Scissor: only [4, 8) x [2, 5) is touched. */
    ASSERT_TRUE(r_begin(&f, c, MD_GPU_LOAD_CLEAR, 0, 0, 0, 1, NULL, 0));
    md_gpu_rect_t sc = {4, 2, 4, 3};
    md_gpu_set_scissor(f.gfx, &sc);
    rargs_t a = {0};
    a.pos = vb;
    a.color = f4(1, 0, 0, 1);
    EXPECT_TRUE(MD_GPU_DRAW(f.gfx, p, 3, 1, a));
    ASSERT_TRUE(md_gpu_render_end(f.gfx));
    ASSERT_TRUE(r_read(&f, c, 0, 0, px, sizeof(px)));
    uint32_t inside = 0, outside = 0;
    for (int y = 0; y < H; ++y) for (int x = 0; x < W; ++x) {
        const bool in = x >= 4 && x < 8 && y >= 2 && y < 5;
        if (in && px[y * W + x] == red) ++inside;
        if (!in && px[y * W + x] == black) ++outside;
    }
    EXPECT_EQ(12u, inside);
    EXPECT_EQ((uint32_t)(W * H - 12), outside);

    /* Viewport: the right half, top-left triangle lands in its top-left.
       A scissor far outside the target is clamped, not an error. */
    ASSERT_TRUE(r_begin(&f, c, MD_GPU_LOAD_CLEAR, 0, 0, 0, 1, NULL, 0));
    md_gpu_rect_t big = {0, 0, 1000, 1000};
    md_gpu_set_scissor(f.gfx, &big);
    md_gpu_viewport_t vp = {8, 0, 8, 16, 0, 1};
    md_gpu_set_viewport(f.gfx, &vp);
    a.pos = vb + 3 * sizeof(md_gpu_float4);
    EXPECT_TRUE(MD_GPU_DRAW(f.gfx, p, 3, 1, a));
    ASSERT_TRUE(md_gpu_render_end(f.gfx));
    ASSERT_TRUE(r_read(&f, c, 0, 0, px, sizeof(px)));
    EXPECT_EQ(red,   px[1 * W + 9]);      /* top-left of the viewport  */
    EXPECT_EQ(black, px[1 * W + 1]);      /* left of the viewport      */
    EXPECT_EQ(black, px[14 * W + 9]);     /* bottom of the viewport    */

    md_gpu_free(f.gfx, vb);
    md_gpu_pipeline_destroy(p);
    md_gpu_texture_destroy(c);
    r_close(&f);
}

/* =========================================================================
   Indexed, strips with restart, points and lines
   ========================================================================= */

UTEST(gpu_render, indexed_u16_u32_and_strip_restart) {
    rfix_t f;
    if (!r_open(&f)) UTEST_SKIP(r_skip_reason());

    md_gpu_texture_t c = r_target(&f, MD_GPU_FORMAT_RGBA8_UNORM, W, H, 0);
    const md_gpu_format_t fmt = MD_GPU_FORMAT_RGBA8_UNORM;
    md_gpu_pipeline_t tris  = r_pipeline(&f, md_shader_gpu_render_test_basic_vs_shader(), md_shader_gpu_render_test_color_fs_shader(),
                                         MD_GPU_TOPOLOGY_TRIANGLES, &fmt, 1, MD_GPU_FORMAT_INVALID);
    md_gpu_pipeline_t strip = r_pipeline(&f, md_shader_gpu_render_test_basic_vs_shader(), md_shader_gpu_render_test_color_fs_shader(),
                                         MD_GPU_TOPOLOGY_TRIANGLE_STRIP, &fmt, 1, MD_GPU_FORMAT_INVALID);
    ASSERT_TRUE(c && tris && strip);

    /* Left half quad (0..3) and right half quad (4..7). */
    const md_gpu_float4 pos[8] = {
        {-1, -1, 0, 1}, {0, -1, 0, 1}, {-1, 1, 0, 1}, {0, 1, 0, 1},
        { 0, -1, 0, 1}, {1, -1, 0, 1}, { 0, 1, 0, 1}, {1, 1, 0, 1},
    };
    const uint16_t i16[6] = {0, 1, 2, 2, 1, 3};
    const uint32_t i32[6] = {4, 5, 6, 6, 5, 7};
    const uint16_t s16[9] = {0, 1, 2, 3, 0xFFFF, 4, 5, 6, 7};
    md_gpu_addr_t vb = r_buffer(&f, pos, sizeof(pos));
    md_gpu_addr_t ib16 = r_buffer(&f, i16, sizeof(i16));
    md_gpu_addr_t ib32 = r_buffer(&f, i32, sizeof(i32));
    md_gpu_addr_t sb16 = r_buffer(&f, s16, sizeof(s16));
    static uint32_t px[W * H];
    const uint32_t red = rgba8(255, 0, 0, 255), blue = rgba8(0, 0, 255, 255);

    ASSERT_TRUE(r_begin(&f, c, MD_GPU_LOAD_CLEAR, 0, 0, 0, 1, NULL, 0));
    rargs_t a = {0};
    a.pos = vb;
    a.color = f4(1, 0, 0, 1);
    EXPECT_TRUE(MD_GPU_DRAW_INDEXED(f.gfx, tris, ib16, MD_GPU_INDEX_U16, 6, 1, a));
    a.color = f4(0, 0, 1, 1);
    EXPECT_TRUE(MD_GPU_DRAW_INDEXED(f.gfx, tris, ib32, MD_GPU_INDEX_U32, 6, 1, a));
    ASSERT_TRUE(md_gpu_render_end(f.gfx));
    ASSERT_TRUE(r_read(&f, c, 0, 0, px, sizeof(px)));
    for (int y = 0; y < H; ++y) for (int x = 0; x < W; ++x) EXPECT_EQ(x < 8 ? red : blue, px[y * W + x]);

    /* One strip, restarted between the quads: both covered, and nothing
       joins them (a joined strip would still cover, so check the count of
       covered pixels via a different colour on a cleared target). */
    ASSERT_TRUE(r_begin(&f, c, MD_GPU_LOAD_CLEAR, 0, 0, 0, 1, NULL, 0));
    a.color = f4(0, 1, 0, 1);
    EXPECT_TRUE(MD_GPU_DRAW_INDEXED(f.gfx, strip, sb16, MD_GPU_INDEX_U16, 9, 1, a));
    ASSERT_TRUE(md_gpu_render_end(f.gfx));
    ASSERT_TRUE(r_read(&f, c, 0, 0, px, sizeof(px)));
    for (int i = 0; i < W * H; ++i) EXPECT_EQ(rgba8(0, 255, 0, 255), px[i]);

    /* Misaligned index address. */
    EXPECT_TRUE(r_begin(&f, c, MD_GPU_LOAD_DONT_CARE, 0, 0, 0, 0, NULL, 0));
    EXPECT_FALSE(MD_GPU_DRAW_INDEXED(f.gfx, tris, ib32 + 2, MD_GPU_INDEX_U32, 3, 1, a));
    EXPECT_TRUE(md_gpu_render_end(f.gfx));

    md_gpu_free(f.gfx, vb); md_gpu_free(f.gfx, ib16); md_gpu_free(f.gfx, ib32); md_gpu_free(f.gfx, sb16);
    md_gpu_pipeline_destroy(tris);
    md_gpu_pipeline_destroy(strip);
    md_gpu_texture_destroy(c);
    r_close(&f);
}

UTEST(gpu_render, points_and_lines) {
    rfix_t f;
    if (!r_open(&f)) UTEST_SKIP(r_skip_reason());

    md_gpu_texture_t c = r_target(&f, MD_GPU_FORMAT_RGBA8_UNORM, W, H, 0);
    const md_gpu_format_t fmt = MD_GPU_FORMAT_RGBA8_UNORM;
    md_gpu_pipeline_t pts = r_pipeline(&f, md_shader_gpu_render_test_point_vs_shader(), md_shader_gpu_render_test_color_fs_shader(),
                                       MD_GPU_TOPOLOGY_POINTS, &fmt, 1, MD_GPU_FORMAT_INVALID);
    md_gpu_pipeline_t lines = r_pipeline(&f, md_shader_gpu_render_test_basic_vs_shader(), md_shader_gpu_render_test_color_fs_shader(),
                                         MD_GPU_TOPOLOGY_LINES, &fmt, 1, MD_GPU_FORMAT_INVALID);
    ASSERT_TRUE(c && pts && lines);

    /* Points at the centres of pixels (2, 3) and (12, 9); a horizontal line
       through the centres of row 13, from x = 0 to x = 16. */
    const md_gpu_float4 pos[4] = {
        {-1 + 2.5f / 8, 1 - 3.5f / 8, 0, 1}, {-1 + 12.5f / 8, 1 - 9.5f / 8, 0, 1},
        {-1, 1 - 13.5f / 8, 0, 1}, {1, 1 - 13.5f / 8, 0, 1},
    };
    md_gpu_addr_t vb = r_buffer(&f, pos, sizeof(pos));
    ASSERT_TRUE(r_begin(&f, c, MD_GPU_LOAD_CLEAR, 0, 0, 0, 1, NULL, 0));
    rargs_t a = {0};
    a.pos = vb;
    a.color = f4(1, 1, 0, 1);
    EXPECT_TRUE(MD_GPU_DRAW(f.gfx, pts, 2, 1, a));
    a.pos = vb + 2 * sizeof(md_gpu_float4);
    a.color = f4(0, 1, 1, 1);
    EXPECT_TRUE(MD_GPU_DRAW(f.gfx, lines, 2, 1, a));
    ASSERT_TRUE(md_gpu_render_end(f.gfx));

    static uint32_t px[W * H];
    ASSERT_TRUE(r_read(&f, c, 0, 0, px, sizeof(px)));
    const uint32_t yellow = rgba8(255, 255, 0, 255), cyan = rgba8(0, 255, 255, 255);
    uint32_t n_yellow = 0, n_cyan_row = 0;
    for (int i = 0; i < W * H; ++i) n_yellow += px[i] == yellow;
    for (int x = 0; x < W; ++x) n_cyan_row += px[13 * W + x] == cyan;
    EXPECT_EQ(yellow, px[3 * W + 2]);
    EXPECT_EQ(yellow, px[9 * W + 12]);
    EXPECT_EQ(2u, n_yellow);
    EXPECT_GE(n_cyan_row, 15u);

    md_gpu_free(f.gfx, vb);
    md_gpu_pipeline_destroy(pts);
    md_gpu_pipeline_destroy(lines);
    md_gpu_texture_destroy(c);
    r_close(&f);
}

/* =========================================================================
   Indirect
   ========================================================================= */

UTEST(gpu_render, multi_draw_indirect_with_first_instance_records) {
    rfix_t f;
    if (!r_open(&f)) UTEST_SKIP(r_skip_reason());

    md_gpu_texture_t u = r_target(&f, MD_GPU_FORMAT_R32_UINT, W, H, 0);
    const md_gpu_format_t fmt = MD_GPU_FORMAT_R32_UINT;
    md_gpu_pipeline_t p = r_pipeline(&f, md_shader_gpu_render_test_indirect_vs_shader(), md_shader_gpu_render_test_id_fs_shader(),
                                     MD_GPU_TOPOLOGY_TRIANGLES, &fmt, 1, MD_GPU_FORMAT_INVALID);
    ASSERT_TRUE(u && p);

    /* Three triangles, each covering one quadrant's top-left half: TL, TR, BL. */
    const md_gpu_float4 pos[9] = {
        {-1, 1, 0, 1}, {-1, 0, 0, 1}, {0, 1, 0, 1},
        { 0, 1, 0, 1}, { 0, 0, 0, 1}, {1, 1, 0, 1},
        {-1, 0, 0, 1}, {-1,-1, 0, 1}, {0, 0, 0, 1},
    };
    rdraw_t rec[3];
    memset(rec, 0, sizeof(rec));
    for (int i = 0; i < 3; ++i) rec[i].color = f4(1, 1, 1, 1);
    const md_gpu_draw_cmd_t cmds[3] = {
        {3, 1, 0, 0},       /* TL, record 0                      */
        {3, 2, 3, 1},       /* TR, record 1, two instances        */
        {3, 0, 6, 2},       /* BL, culled: zero instances         */
    };
    md_gpu_addr_t vb = r_buffer(&f, pos, sizeof(pos));
    md_gpu_addr_t rb = r_buffer(&f, rec, sizeof(rec));
    md_gpu_addr_t cb = r_buffer(&f, cmds, sizeof(cmds));

    md_gpu_render_desc_t rd = {0};
    rd.color_count = 1;
    rd.color[0].texture = u; rd.color[0].load = MD_GPU_LOAD_CLEAR; rd.color[0].clear.u32[0] = 0xFFFFFFFFu;
    ASSERT_TRUE(md_gpu_render_begin(f.gfx, &rd));
    rargs_t a = {0};
    a.pos   = vb;
    a.draws = rb;
    a.flags = RF_RECORDS;
    a.id    = 1000;
    EXPECT_TRUE(md_gpu_draw_indirect(f.gfx, p, cb, 3, &a, sizeof(a)));
    ASSERT_TRUE(md_gpu_render_end(f.gfx));

    static uint32_t ids[W * H];
    ASSERT_TRUE(r_read(&f, u, 0, 0, ids, sizeof(ids)));
    EXPECT_EQ(1000u + 0u,  ids[1 * W + 1]);          /* record 0, instance 0 */
    EXPECT_EQ(1000u + 17u, ids[1 * W + 9]);          /* record 1, instance 1 */
    EXPECT_EQ(0xFFFFFFFFu, ids[9 * W + 1]);          /* culled               */

    /* Indexed indirect: vertex_offset and first_index both apply. */
    const uint32_t idx[6] = {99, 99, 99, 0, 1, 2};
    const md_gpu_draw_indexed_cmd_t icmds[1] = {{3, 1, 3, 6, 2}};   /* indices 3..5, +6 -> BL, record 2 */
    md_gpu_addr_t ib = r_buffer(&f, idx, sizeof(idx));
    md_gpu_addr_t icb = r_buffer(&f, icmds, sizeof(icmds));
    ASSERT_TRUE(md_gpu_render_begin(f.gfx, &rd));
    EXPECT_TRUE(md_gpu_draw_indexed_indirect(f.gfx, p, ib, MD_GPU_INDEX_U32, icb, 1, &a, sizeof(a)));
    ASSERT_TRUE(md_gpu_render_end(f.gfx));
    ASSERT_TRUE(r_read(&f, u, 0, 0, ids, sizeof(ids)));
    EXPECT_EQ(1000u + 32u, ids[9 * W + 1]);          /* record 2, instance 0 */
    EXPECT_EQ(0xFFFFFFFFu, ids[1 * W + 1]);

    md_gpu_free(f.gfx, vb); md_gpu_free(f.gfx, rb); md_gpu_free(f.gfx, cb);
    md_gpu_free(f.gfx, ib); md_gpu_free(f.gfx, icb);
    md_gpu_pipeline_destroy(p);
    md_gpu_texture_destroy(u);
    r_close(&f);
}

/* =========================================================================
   Attachments at mips and layers; ordering between passes and kernels
   ========================================================================= */

UTEST(gpu_render, mip_and_layer_attachments) {
    rfix_t f;
    if (!r_open(&f)) UTEST_SKIP(r_skip_reason());

    md_gpu_texture_desc_t td = {0};
    td.type   = MD_GPU_TEX_2D_ARRAY;
    td.format = MD_GPU_FORMAT_RGBA8_UNORM;
    td.usage  = MD_GPU_TEX_RENDER_TARGET | MD_GPU_TEX_SAMPLED;
    td.width  = W; td.height = H; td.depth_or_layers = 3; td.mip_levels = 2;
    md_gpu_texture_t t = md_gpu_texture_create(f.gfx, &td);
    ASSERT_TRUE(t != NULL);

    for (uint32_t mip = 0; mip < 2; ++mip) for (uint32_t layer = 0; layer < 3; ++layer) {
        md_gpu_render_desc_t rd = {0};
        rd.color_count = 1;
        rd.color[0].texture = t;
        rd.color[0].mip     = mip;
        rd.color[0].layer   = layer;
        rd.color[0].load    = MD_GPU_LOAD_CLEAR;
        rd.color[0].clear.f32[0] = (float)mip;
        rd.color[0].clear.f32[1] = (float)layer / 2.0f;
        rd.color[0].clear.f32[3] = 1.0f;
        ASSERT_TRUE(md_gpu_render_begin(f.gfx, &rd));
        ASSERT_TRUE(md_gpu_render_end(f.gfx));
    }
    static uint32_t px[W * H];
    ASSERT_TRUE(r_read(&f, t, 1, 2, px, (W / 2) * (H / 2) * 4));
    EXPECT_EQ(rgba8(255, 255, 0, 255), px[0]);
    EXPECT_EQ(rgba8(255, 255, 0, 255), px[(W / 2) * (H / 2) - 1]);
    ASSERT_TRUE(r_read(&f, t, 0, 1, px, sizeof(px)));
    EXPECT_TRUE(near8(px[0], rgba8(0, 128, 0, 255), 1));

    /* Out-of-range layer and mip are rejected. */
    md_gpu_render_desc_t bad = {0};
    bad.color_count = 1;
    bad.color[0].texture = t;
    bad.color[0].layer   = 3;
    EXPECT_FALSE(md_gpu_render_begin(f.gfx, &bad));
    bad.color[0].layer = 0;
    bad.color[0].mip   = 2;
    EXPECT_FALSE(md_gpu_render_begin(f.gfx, &bad));

    md_gpu_texture_destroy(t);
    r_close(&f);
}

UTEST(gpu_render, attachment_sizes_must_agree) {
    rfix_t f;
    if (!r_open(&f)) UTEST_SKIP(r_skip_reason());
    md_gpu_texture_t a = r_target(&f, MD_GPU_FORMAT_RGBA8_UNORM, W, H, 0);
    md_gpu_texture_t d = r_target(&f, MD_GPU_FORMAT_D32_FLOAT, W / 2, H, 0);
    ASSERT_TRUE(a && d);
    EXPECT_FALSE(r_begin(&f, a, MD_GPU_LOAD_CLEAR, 0, 0, 0, 0, d, 1.0f));
    /* A depth format as a colour attachment, and the reverse. */
    EXPECT_FALSE(r_begin(&f, d, MD_GPU_LOAD_CLEAR, 0, 0, 0, 0, NULL, 1.0f));
    EXPECT_FALSE(r_begin(&f, NULL, MD_GPU_LOAD_CLEAR, 0, 0, 0, 0, a, 1.0f));
    md_gpu_texture_destroy(a);
    md_gpu_texture_destroy(d);
    r_close(&f);
}

/* Pass 1 renders A; pass 2 samples A into B. In EXPLICIT mode the caller
   states the dependency; in IMPLICIT mode the stream does. */
static void render_then_sample(rfix_t* f, bool explicit_mode, int* utest_result) {
    md_gpu_texture_t ta = r_target(f, MD_GPU_FORMAT_RGBA8_UNORM, W, H, MD_GPU_TEX_SAMPLED);
    md_gpu_texture_t tb = r_target(f, MD_GPU_FORMAT_RGBA8_UNORM, W, H, 0);
    const md_gpu_format_t fmt = MD_GPU_FORMAT_RGBA8_UNORM;
    md_gpu_pipeline_t p = r_pipeline(f, md_shader_gpu_render_test_basic_vs_shader(), md_shader_gpu_render_test_color_fs_shader(),
                                     MD_GPU_TOPOLOGY_TRIANGLES, &fmt, 1, MD_GPU_FORMAT_INVALID);
    md_gpu_addr_t vb = r_buffer(f, full_tri, sizeof(full_tri));
    md_gpu_sampler_desc_t sd = {0};
    md_gpu_sampler_t smp = md_gpu_sampler(f->dev, &sd);
    ASSERT_TRUE(ta && tb && p && vb && smp.handle);

    /* Entering EXPLICIT orders the region after the setup work above. */
    if (explicit_mode) md_gpu_stream_set_ordering(f->gfx, MD_GPU_ORDER_EXPLICIT);
    rargs_t a = {0};
    a.pos = vb;
    for (int frame = 0; frame < 3; ++frame) {
        const float g = 0.25f * (float)(frame + 1);
        ASSERT_TRUE(r_begin(f, ta, MD_GPU_LOAD_DONT_CARE, 0, 0, 0, 0, NULL, 0));
        a.color = f4(0, g, 0, 1);
        a.flags = 0;
        EXPECT_TRUE(MD_GPU_DRAW(f->gfx, p, 3, 1, a));
        ASSERT_TRUE(md_gpu_render_end(f->gfx));
        if (explicit_mode) md_gpu_barrier(f->gfx, MD_GPU_STAGE_ATTACHMENT, MD_GPU_STAGE_FRAGMENT);

        ASSERT_TRUE(r_begin(f, tb, MD_GPU_LOAD_DONT_CARE, 0, 0, 0, 0, NULL, 0));
        a.color = f4(1, 1, 1, 1);
        a.flags = RF_SAMPLE;
        a.tex   = md_gpu_texture_sampled(ta);
        a.smp   = smp;
        EXPECT_TRUE(MD_GPU_DRAW(f->gfx, p, 3, 1, a));
        ASSERT_TRUE(md_gpu_render_end(f->gfx));
        if (explicit_mode) md_gpu_barrier(f->gfx, MD_GPU_STAGE_ATTACHMENT | MD_GPU_STAGE_FRAGMENT,
                                          MD_GPU_STAGE_TRANSFER | MD_GPU_STAGE_ATTACHMENT);

        static uint32_t px[W * H];
        ASSERT_TRUE(r_read(f, tb, 0, 0, px, sizeof(px)));
        const uint8_t gg = (uint8_t)(g * 255.0f + 0.5f);
        EXPECT_TRUE(near8(px[0], rgba8(0, gg, 0, 255), 1));
        EXPECT_TRUE(near8(px[W * H - 1], rgba8(0, gg, 0, 255), 1));
        /* The readback copy is the last reader of `tb` before the next
           frame's pass writes it. */
        if (explicit_mode) md_gpu_barrier(f->gfx, MD_GPU_STAGE_TRANSFER, MD_GPU_STAGE_ATTACHMENT | MD_GPU_STAGE_TRANSFER);
    }
    md_gpu_stream_set_ordering(f->gfx, MD_GPU_ORDER_IMPLICIT);
    md_gpu_free(f->gfx, vb);
    md_gpu_pipeline_destroy(p);
    md_gpu_texture_destroy(ta);
    md_gpu_texture_destroy(tb);
}

UTEST(gpu_render, pass_to_pass_implicit) {
    rfix_t f;
    if (!r_open(&f)) UTEST_SKIP(r_skip_reason());
    render_then_sample(&f, false, utest_result);
    r_close(&f);
}

UTEST(gpu_render, pass_to_pass_explicit_barriers) {
    rfix_t f;
    if (!r_open(&f)) UTEST_SKIP(r_skip_reason());
    render_then_sample(&f, true, utest_result);
    r_close(&f);
}

UTEST(gpu_render, kernel_feeds_draw_and_reads_the_result) {
    rfix_t f;
    if (!r_open(&f)) UTEST_SKIP(r_skip_reason());

    md_gpu_kernel_desc_t kd = md_shader_gpu_render_test_gen_pos_kernel();
    md_gpu_kernel_t gen = md_gpu_kernel_create(f.dev, &kd);
    kd = md_shader_gpu_render_test_read_tex_kernel();
    md_gpu_kernel_t read = md_gpu_kernel_create(f.dev, &kd);
    md_gpu_texture_t c = r_target(&f, MD_GPU_FORMAT_RGBA32_FLOAT, W, H, MD_GPU_TEX_SAMPLED);
    const md_gpu_format_t fmt = MD_GPU_FORMAT_RGBA32_FLOAT;
    md_gpu_pipeline_t p = r_pipeline(&f, md_shader_gpu_render_test_basic_vs_shader(), md_shader_gpu_render_test_color_fs_shader(),
                                     MD_GPU_TOPOLOGY_TRIANGLES, &fmt, 1, MD_GPU_FORMAT_INVALID);
    ASSERT_TRUE(gen && read && c && p);

    /* The kernel scales a half-size triangle up to the full-screen one. */
    md_gpu_float4 half[3];
    for (int v = 0; v < 3; ++v) half[v] = f4(full_tri[v].x * 0.5f, full_tri[v].y * 0.5f, 0.25f, 0.5f);
    md_gpu_addr_t src = r_buffer(&f, half, sizeof(half));
    md_gpu_addr_t pos = r_buffer(&f, NULL, sizeof(half));
    md_gpu_addr_t out = r_buffer(&f, NULL, W * H * sizeof(md_gpu_float4));

    for (int round = 0; round < 2; ++round) {
        rargs_t a = {0};
        a.color  = f4(2, 2, 2, 2);
        a.pos    = pos;
        a.colors = src;
        a.width  = 3;
        EXPECT_TRUE(md_gpu_launch(f.gfx, gen, md_gpu_grid(1, 1, 1), &a, sizeof(a)));

        ASSERT_TRUE(r_begin(&f, c, MD_GPU_LOAD_CLEAR, 0, 0, 0, 0, NULL, 0));
        a.color = f4(0.5f, 0.25f, (float)round, 1);
        EXPECT_TRUE(MD_GPU_DRAW(f.gfx, p, 3, 1, a));
        ASSERT_TRUE(md_gpu_render_end(f.gfx));

        rargs_t r = {0};
        r.tex    = md_gpu_texture_sampled(c);
        r.colors = out;
        r.width  = W;
        r.height = H;
        EXPECT_TRUE(md_gpu_launch(f.gfx, read, md_gpu_grid(W / 8, H / 8, 1), &r, sizeof(r)));

        static md_gpu_float4 host[W * H];
        md_gpu_mem_t rb = md_gpu_malloc(f.gfx, MD_GPU_MEM_HOST_READ, sizeof(host));
        md_gpu_copy(f.gfx, rb.gpu, out, sizeof(host));
        md_gpu_stream_sync(f.gfx);
        memcpy(host, rb.cpu, sizeof(host));
        md_gpu_free(f.gfx, rb.gpu);
        EXPECT_EQ(0.5f, host[0].x);
        EXPECT_EQ(0.25f, host[W * H - 1].y);
        EXPECT_EQ((float)round, host[W * 3 + 11].z);
    }

    md_gpu_free(f.gfx, src); md_gpu_free(f.gfx, pos); md_gpu_free(f.gfx, out);
    md_gpu_kernel_destroy(gen);
    md_gpu_kernel_destroy(read);
    md_gpu_pipeline_destroy(p);
    md_gpu_texture_destroy(c);
    r_close(&f);
}

/* Rendered on the graphics stream, read on a compute stream that joins it. */
UTEST(gpu_render, compute_stream_reads_a_rendered_texture) {
    rfix_t f;
    if (!r_open(&f)) UTEST_SKIP(r_skip_reason());
    md_gpu_stream_t comp = md_gpu_stream_create(f.dev, MD_GPU_STREAM_COMPUTE, "reader");
    md_gpu_kernel_desc_t kd = md_shader_gpu_render_test_read_tex_kernel();
    md_gpu_kernel_t read = md_gpu_kernel_create(f.dev, &kd);
    md_gpu_texture_t c = r_target(&f, MD_GPU_FORMAT_RGBA32_FLOAT, W, H, MD_GPU_TEX_SAMPLED);
    ASSERT_TRUE(comp && read && c);

    md_gpu_render_desc_t rd = {0};
    rd.color_count = 1;
    rd.color[0].texture = c;
    rd.color[0].load    = MD_GPU_LOAD_CLEAR;
    rd.color[0].clear.f32[0] = 3.0f; rd.color[0].clear.f32[3] = 7.0f;
    ASSERT_TRUE(md_gpu_render_begin(f.gfx, &rd));
    ASSERT_TRUE(md_gpu_render_end(f.gfx));
    md_gpu_stream_wait(comp, md_gpu_stream_record(f.gfx));

    md_gpu_addr_t out = md_gpu_malloc(comp, MD_GPU_MEM_DEVICE, W * H * sizeof(md_gpu_float4)).gpu;
    rargs_t r = {0};
    r.tex = md_gpu_texture_sampled(c);
    r.colors = out;
    r.width = W; r.height = H;
    EXPECT_TRUE(md_gpu_launch(comp, read, md_gpu_grid(W / 8, H / 8, 1), &r, sizeof(r)));
    md_gpu_mem_t rb = md_gpu_malloc(comp, MD_GPU_MEM_HOST_READ, W * H * sizeof(md_gpu_float4));
    md_gpu_copy(comp, rb.gpu, out, W * H * sizeof(md_gpu_float4));
    md_gpu_stream_sync(comp);
    const md_gpu_float4* host = (const md_gpu_float4*)rb.cpu;
    EXPECT_EQ(3.0f, host[5].x);
    EXPECT_EQ(7.0f, host[W * H - 1].w);

    md_gpu_free(comp, out); md_gpu_free(comp, rb.gpu);
    md_gpu_kernel_destroy(read);
    md_gpu_texture_destroy(c);
    md_gpu_stream_destroy(comp);
    r_close(&f);
}

/* =========================================================================
   Texture to texture copies
   ========================================================================= */

UTEST(gpu_render, copy_texture_subregion) {
    rfix_t f;
    if (!r_open(&f)) UTEST_SKIP(r_skip_reason());

    md_gpu_texture_desc_t td = {0};
    td.type = MD_GPU_TEX_2D; td.format = MD_GPU_FORMAT_R32_UINT; td.usage = MD_GPU_TEX_SAMPLED;
    td.width = 8; td.height = 8;
    md_gpu_texture_t src = md_gpu_texture_create(f.gfx, &td);
    td.width = 16; td.height = 16;
    md_gpu_texture_t dst = md_gpu_texture_create(f.gfx, &td);
    td.format = MD_GPU_FORMAT_R32_FLOAT;
    md_gpu_texture_t other = md_gpu_texture_create(f.gfx, &td);
    ASSERT_TRUE(src && dst && other);

    uint32_t s[64], z[256];
    for (int i = 0; i < 64; ++i) s[i] = (uint32_t)i;
    memset(z, 0, sizeof(z));
    ASSERT_TRUE(md_gpu_upload_texture(f.gfx, src, NULL, s, sizeof(s)));
    ASSERT_TRUE(md_gpu_upload_texture(f.gfx, dst, NULL, z, sizeof(z)));

    md_gpu_tex_region_t sr = {0}, dr = {0};
    sr.offset[0] = 2; sr.offset[1] = 1; sr.extent[0] = 4; sr.extent[1] = 3;
    dr.offset[0] = 10; dr.offset[1] = 12;
    ASSERT_TRUE(md_gpu_copy_texture(f.gfx, dst, &dr, src, &sr));

    static uint32_t out[256];
    ASSERT_TRUE(r_read(&f, dst, 0, 0, out, sizeof(out)));
    for (int y = 0; y < 16; ++y) for (int x = 0; x < 16; ++x) {
        const bool in = x >= 10 && x < 14 && y >= 12 && y < 15;
        const uint32_t expect = in ? (uint32_t)((y - 12 + 1) * 8 + (x - 10 + 2)) : 0u;
        EXPECT_EQ(expect, out[y * 16 + x]);
    }

    /* Format mismatch and a destination overrun fail. */
    EXPECT_FALSE(md_gpu_copy_texture(f.gfx, other, NULL, src, NULL));
    dr.offset[0] = 14;
    EXPECT_FALSE(md_gpu_copy_texture(f.gfx, dst, &dr, src, &sr));

    md_gpu_texture_destroy(src);
    md_gpu_texture_destroy(dst);
    md_gpu_texture_destroy(other);
    r_close(&f);
}

/* =========================================================================
   Rules and errors
   ========================================================================= */

static void r_count_fn(void* user) { ++*(int*)user; }

UTEST(gpu_render, inside_a_pass_only_draws_and_state_are_allowed) {
    rfix_t f;
    if (!r_open(&f)) UTEST_SKIP(r_skip_reason());

    md_gpu_texture_t c = r_target(&f, MD_GPU_FORMAT_RGBA8_UNORM, W, H, 0);
    md_gpu_format_t fmt = MD_GPU_FORMAT_RGBA8_UNORM;
    md_gpu_pipeline_t p = r_pipeline(&f, md_shader_gpu_render_test_basic_vs_shader(), md_shader_gpu_render_test_color_fs_shader(),
                                     MD_GPU_TOPOLOGY_TRIANGLES, &fmt, 1, MD_GPU_FORMAT_INVALID);
    fmt = MD_GPU_FORMAT_RGBA16_FLOAT;
    md_gpu_pipeline_t p16 = r_pipeline(&f, md_shader_gpu_render_test_basic_vs_shader(), md_shader_gpu_render_test_color_fs_shader(),
                                       MD_GPU_TOPOLOGY_TRIANGLES, &fmt, 1, MD_GPU_FORMAT_INVALID);
    md_gpu_pipeline_t pd = r_pipeline(&f, md_shader_gpu_render_test_basic_vs_shader(), md_shader_gpu_render_test_color_fs_shader(),
                                      MD_GPU_TOPOLOGY_TRIANGLES, &(md_gpu_format_t){MD_GPU_FORMAT_RGBA8_UNORM}, 1, MD_GPU_FORMAT_D32_FLOAT);
    md_gpu_addr_t vb = r_buffer(&f, full_tri, sizeof(full_tri));
    md_gpu_addr_t buf = r_buffer(&f, NULL, 256);
    md_gpu_stream_t comp = md_gpu_stream_default(f.dev, MD_GPU_STREAM_COMPUTE);
    ASSERT_TRUE(c && p && p16 && pd && vb && buf);

    rargs_t a = {0};
    a.pos = vb;
    EXPECT_FALSE(MD_GPU_DRAW(f.gfx, p, 3, 1, a));                 /* outside a pass */
    md_gpu_set_draw_state(f.gfx, NULL);                           /* error, harmless */

    ASSERT_TRUE(r_begin(&f, c, MD_GPU_LOAD_CLEAR, 0, 0, 0, 1, NULL, 0));
    EXPECT_FALSE(r_begin(&f, c, MD_GPU_LOAD_CLEAR, 0, 0, 0, 1, NULL, 0));   /* nested */
    uint32_t word = 5;
    EXPECT_FALSE(md_gpu_upload(f.gfx, buf, &word, 4));
    EXPECT_FALSE(md_gpu_memset(f.gfx, buf, 0, 16));
    EXPECT_FALSE(md_gpu_copy(f.gfx, buf, buf + 128, 16));
    EXPECT_TRUE(md_gpu_upload_begin(f.gfx, buf, 4) == NULL);
    md_gpu_texture_desc_t td = {0};
    td.type = MD_GPU_TEX_2D; td.format = MD_GPU_FORMAT_R32_FLOAT; td.usage = MD_GPU_TEX_SAMPLED; td.width = td.height = 4;
    EXPECT_TRUE(md_gpu_texture_create(f.gfx, &td) == NULL);
    EXPECT_FALSE(md_gpu_sync_is_valid(md_gpu_stream_record(f.gfx)));
    md_gpu_stream_wait(f.gfx, md_gpu_stream_record(comp));
    int fired = 0;
    EXPECT_FALSE(md_gpu_launch_host_fn(f.gfx, r_count_fn, &fired));
    /* Memory and temp scopes record nothing and stay legal. */
    md_gpu_mem_t m = md_gpu_malloc(f.gfx, MD_GPU_MEM_DEVICE, 64);
    EXPECT_TRUE(m.gpu != 0);
    md_gpu_temp_t scope = md_gpu_temp_begin(f.gfx);
    EXPECT_TRUE(md_gpu_temp_alloc(f.gfx, MD_GPU_MEM_HOST_WRITE, 64).cpu != NULL);
    md_gpu_temp_end(f.gfx, scope);
    /* Format mismatches: colour format, and a depth format the pass lacks. */
    EXPECT_FALSE(MD_GPU_DRAW(f.gfx, p16, 3, 1, a));
    EXPECT_FALSE(MD_GPU_DRAW(f.gfx, pd, 3, 1, a));
    /* Argument size mismatch. */
    EXPECT_FALSE(md_gpu_draw(f.gfx, p, 3, 1, &a, sizeof(a) - 16));
    /* The right draw still works after all that. */
    a.color = f4(0, 0, 1, 1);
    EXPECT_TRUE(MD_GPU_DRAW(f.gfx, p, 3, 1, a));
    ASSERT_TRUE(md_gpu_render_end(f.gfx));
    EXPECT_FALSE(md_gpu_render_end(f.gfx));                        /* nothing open */
    md_gpu_free(f.gfx, m.gpu);
    md_gpu_device_poll(f.dev);
    EXPECT_EQ(0, fired);

    static uint32_t px[W * H];
    ASSERT_TRUE(r_read(&f, c, 0, 0, px, sizeof(px)));
    EXPECT_EQ(rgba8(0, 0, 255, 255), px[W * H / 2]);

    /* Passes need a GRAPHICS stream. */
    md_gpu_render_desc_t rd = {0};
    rd.color_count = 1;
    rd.color[0].texture = c;
    EXPECT_FALSE(md_gpu_render_begin(comp, &rd));

    md_gpu_free(f.gfx, vb);
    md_gpu_free(f.gfx, buf);
    md_gpu_pipeline_destroy(p);
    md_gpu_pipeline_destroy(p16);
    md_gpu_pipeline_destroy(pd);
    md_gpu_texture_destroy(c);
    r_close(&f);
}

UTEST(gpu_render, pipeline_and_target_validation) {
    rfix_t f;
    if (!r_open(&f)) UTEST_SKIP(r_skip_reason());

    /* Blending an integer target. */
    md_gpu_pipeline_desc_t pd = {0};
    pd.vertex   = md_shader_gpu_render_test_basic_vs_shader();
    pd.fragment = md_shader_gpu_render_test_id_fs_shader();
    pd.color_count = 1;
    pd.color[0].format = MD_GPU_FORMAT_R32_UINT;
    pd.color[0].blend.enable = true;
    EXPECT_TRUE(md_gpu_pipeline_create(f.dev, &pd) == NULL);
    /* A depth format as a colour target. */
    pd.color[0].format = MD_GPU_FORMAT_D32_FLOAT;
    pd.color[0].blend.enable = false;
    EXPECT_TRUE(md_gpu_pipeline_create(f.dev, &pd) == NULL);
    /* Colour targets without a fragment shader. */
    pd.color[0].format = MD_GPU_FORMAT_RGBA8_UNORM;
    memset(&pd.fragment, 0, sizeof(pd.fragment));
    EXPECT_TRUE(md_gpu_pipeline_create(f.dev, &pd) == NULL);
    /* Stages that disagree on the argument struct. */
    pd.fragment = md_shader_gpu_imgui_test_imgui_fs_shader();
    EXPECT_TRUE(md_gpu_pipeline_create(f.dev, &pd) == NULL);
    /* A non-depth depth format. */
    pd.fragment = md_shader_gpu_render_test_color_fs_shader();
    pd.depth_format = MD_GPU_FORMAT_R32_FLOAT;
    EXPECT_TRUE(md_gpu_pipeline_create(f.dev, &pd) == NULL);

    /* RENDER_TARGET is 2D and 2D_ARRAY only. */
    md_gpu_texture_desc_t td = {0};
    td.type = MD_GPU_TEX_3D; td.format = MD_GPU_FORMAT_RGBA8_UNORM; td.usage = MD_GPU_TEX_RENDER_TARGET;
    td.width = td.height = td.depth_or_layers = 4;
    EXPECT_TRUE(md_gpu_texture_create(f.gfx, &td) == NULL);

    /* A texture without RENDER_TARGET cannot be attached. */
    td.type = MD_GPU_TEX_2D; td.usage = MD_GPU_TEX_SAMPLED; td.depth_or_layers = 1;
    md_gpu_texture_t s = md_gpu_texture_create(f.gfx, &td);
    ASSERT_TRUE(s != NULL);
    EXPECT_FALSE(r_begin(&f, s, MD_GPU_LOAD_CLEAR, 0, 0, 0, 0, NULL, 0));
    md_gpu_texture_destroy(s);
    r_close(&f);
}

/* =========================================================================
   ImGui: the renderer from the design notes, run for several frames with
   per-frame temp memory and two frames in flight
   ========================================================================= */

typedef struct { float pos[2]; float uv[2]; uint32_t col; } ImDrawVert;
typedef uint16_t ImDrawIdx;
typedef struct { float clip_rect[4]; uint64_t texture_id; uint32_t vtx_offset, idx_offset, elem_count; } ImDrawCmd;
typedef struct { const ImDrawCmd* cmds; int cmd_count; const ImDrawVert* vtx; int vtx_count; const ImDrawIdx* idx; int idx_count; } ImDrawList;
typedef struct { const ImDrawList* lists; int list_count; float display_pos[2], display_size[2], fb_scale[2]; } ImDrawData;

/* Mirrors Args in unittest/shaders/gpu_imgui_test.slang. */
typedef struct {
    md_gpu_float2        scale;
    md_gpu_float2        translate;
    md_gpu_addr_t        verts;
    md_gpu_sampled_tex_t tex;        /* ImTextureID is this handle */
    md_gpu_sampler_t     smp;
} imgui_args_t;

typedef struct {
    md_gpu_pipeline_t pipeline;
    md_gpu_sampler_t  sampler;
} imgui_renderer_t;

static bool imgui_renderer_init(imgui_renderer_t* r, md_gpu_device_t dev, md_gpu_format_t target) {
    md_gpu_pipeline_desc_t pd = {0};
    pd.vertex              = md_shader_gpu_imgui_test_imgui_vs_shader();
    pd.fragment            = md_shader_gpu_imgui_test_imgui_fs_shader();
    pd.topology            = MD_GPU_TOPOLOGY_TRIANGLES;
    pd.color_count         = 1;
    pd.color[0].format     = target;
    pd.color[0].blend      = (md_gpu_blend_t){
        .enable    = true,
        .src_color = MD_GPU_BLEND_SRC_ALPHA, .dst_color = MD_GPU_BLEND_ONE_MINUS_SRC_ALPHA,
        .src_alpha = MD_GPU_BLEND_ONE,       .dst_alpha = MD_GPU_BLEND_ONE_MINUS_SRC_ALPHA,
    };
    pd.label = "imgui";
    r->pipeline = md_gpu_pipeline_create(dev, &pd);

    md_gpu_sampler_desc_t sd = {0};
    sd.min_filter = sd.mag_filter = MD_GPU_FILTER_LINEAR;
    r->sampler = md_gpu_sampler(dev, &sd);
    return r->pipeline != NULL;
}

/* Records into an open pass whose colour target 0 has the renderer's format.
   Geometry is temp memory, in a scope of its own nested in the caller's. */
static void imgui_render(const imgui_renderer_t* r, md_gpu_stream_t s, const ImDrawData* dd) {
    const float fb_w = dd->display_size[0] * dd->fb_scale[0];
    const float fb_h = dd->display_size[1] * dd->fb_scale[1];
    if (fb_w <= 0.0f || fb_h <= 0.0f) return;

    imgui_args_t a;
    memset(&a, 0, sizeof(a));
    a.scale.x     =  2.0f / dd->display_size[0];
    a.scale.y     = -2.0f / dd->display_size[1];          /* clip space is +Y up */
    a.translate.x = -1.0f - dd->display_pos[0] * a.scale.x;
    a.translate.y =  1.0f - dd->display_pos[1] * a.scale.y;
    a.smp         = r->sampler;

    const md_gpu_viewport_t vp = { 0.0f, 0.0f, fb_w, fb_h, 0.0f, 1.0f };
    md_gpu_set_viewport(s, &vp);

    md_gpu_temp_t scope = md_gpu_temp_begin(s);
    for (int l = 0; l < dd->list_count; ++l) {
        const ImDrawList* list = &dd->lists[l];
        const size_t vbytes = (size_t)list->vtx_count * sizeof(ImDrawVert);
        const size_t ibytes = (size_t)list->idx_count * sizeof(ImDrawIdx);
        md_gpu_mem_t v = md_gpu_temp_alloc(s, MD_GPU_MEM_HOST_WRITE, vbytes);
        md_gpu_mem_t i = md_gpu_temp_alloc(s, MD_GPU_MEM_HOST_WRITE, ibytes);
        if (!v.cpu || !i.cpu) break;
        memcpy(v.cpu, list->vtx, vbytes);
        memcpy(i.cpu, list->idx, ibytes);

        for (int c = 0; c < list->cmd_count; ++c) {
            const ImDrawCmd* cmd = &list->cmds[c];
            float x0 = (cmd->clip_rect[0] - dd->display_pos[0]) * dd->fb_scale[0];
            float y0 = (cmd->clip_rect[1] - dd->display_pos[1]) * dd->fb_scale[1];
            float x1 = (cmd->clip_rect[2] - dd->display_pos[0]) * dd->fb_scale[0];
            float y1 = (cmd->clip_rect[3] - dd->display_pos[1]) * dd->fb_scale[1];
            if (x0 < 0.0f) x0 = 0.0f;
            if (y0 < 0.0f) y0 = 0.0f;
            if (x1 > fb_w) x1 = fb_w;
            if (y1 > fb_h) y1 = fb_h;
            if (x1 <= x0 || y1 <= y0) continue;
            const md_gpu_rect_t sc = { (uint32_t)x0, (uint32_t)y0, (uint32_t)(x1 - x0), (uint32_t)(y1 - y0) };
            md_gpu_set_scissor(s, &sc);

            a.verts      = v.gpu + (md_gpu_addr_t)cmd->vtx_offset * sizeof(ImDrawVert);
            a.tex.handle = cmd->texture_id;
            MD_GPU_DRAW_INDEXED(s, r->pipeline,
                                i.gpu + (md_gpu_addr_t)cmd->idx_offset * sizeof(ImDrawIdx), MD_GPU_INDEX_U16,
                                cmd->elem_count, 1, a);
        }
    }
    md_gpu_temp_end(s, scope);
}

/* Two rectangles in a 32 x 32 display: a white one over the left half,
   clipped to its top half, and a half-transparent red one over the right
   half, unclipped. */
static void imgui_test_frame(ImDrawData* dd, ImDrawList* list, ImDrawCmd cmds[2], ImDrawVert vtx[8], ImDrawIdx idx[12], uint64_t tex) {
    const float S = 32.0f;
    const uint32_t white = 0xFFFFFFFFu, red_half = 0x800000FFu;
    const float rects[2][4] = {{0, 0, S / 2, S}, {S / 2, 0, S, S}};
    const uint32_t cols[2] = {white, red_half};
    for (int r = 0; r < 2; ++r) {
        const float x0 = rects[r][0], y0 = rects[r][1], x1 = rects[r][2], y1 = rects[r][3];
        ImDrawVert* v = &vtx[r * 4];
        v[0] = (ImDrawVert){{x0, y0}, {0, 0}, cols[r]};
        v[1] = (ImDrawVert){{x1, y0}, {1, 0}, cols[r]};
        v[2] = (ImDrawVert){{x1, y1}, {1, 1}, cols[r]};
        v[3] = (ImDrawVert){{x0, y1}, {0, 1}, cols[r]};
        const ImDrawIdx q[6] = {0, 1, 2, 0, 2, 3};
        memcpy(&idx[r * 6], q, sizeof(q));
        cmds[r].texture_id = tex;
        cmds[r].vtx_offset = (uint32_t)(r * 4);
        cmds[r].idx_offset = (uint32_t)(r * 6);
        cmds[r].elem_count = 6;
    }
    cmds[0].clip_rect[0] = 0;     cmds[0].clip_rect[1] = 0; cmds[0].clip_rect[2] = S / 2; cmds[0].clip_rect[3] = S / 2;
    cmds[1].clip_rect[0] = -100;  cmds[1].clip_rect[1] = -100; cmds[1].clip_rect[2] = 100; cmds[1].clip_rect[3] = 100;
    list->cmds = cmds; list->cmd_count = 2;
    list->vtx = vtx;   list->vtx_count = 8;
    list->idx = idx;   list->idx_count = 12;
    memset(dd, 0, sizeof(*dd));
    dd->lists = list; dd->list_count = 1;
    dd->display_size[0] = S; dd->display_size[1] = S;
    dd->fb_scale[0] = 1; dd->fb_scale[1] = 1;
}

UTEST(gpu_render, imgui_renderer_frames) {
    rfix_t f;
    if (!r_open(&f)) UTEST_SKIP(r_skip_reason());

    imgui_renderer_t ui;
    ASSERT_TRUE(imgui_renderer_init(&ui, f.dev, MD_GPU_FORMAT_RGBA8_UNORM));
    md_gpu_texture_t target = r_target(&f, MD_GPU_FORMAT_RGBA8_UNORM, 32, 32, 0);
    md_gpu_texture_desc_t td = {0};
    td.type = MD_GPU_TEX_2D; td.format = MD_GPU_FORMAT_RGBA8_UNORM; td.usage = MD_GPU_TEX_SAMPLED;
    td.width = td.height = 2; td.label = "font atlas";
    md_gpu_texture_t font = md_gpu_texture_create(f.gfx, &td);
    ASSERT_TRUE(target && font);
    const uint32_t white4[4] = {0xFFFFFFFFu, 0xFFFFFFFFu, 0xFFFFFFFFu, 0xFFFFFFFFu};
    ASSERT_TRUE(md_gpu_upload_texture(f.gfx, font, NULL, white4, sizeof(white4)));

    ImDrawData dd; ImDrawList list; ImDrawCmd cmds[2]; ImDrawVert vtx[8]; ImDrawIdx idx[12];
    memset(cmds, 0, sizeof(cmds));
    imgui_test_frame(&dd, &list, cmds, vtx, idx, md_gpu_texture_sampled(font).handle);

    md_gpu_sync_t frame_done[2] = {md_gpu_sync_none(), md_gpu_sync_none()};
    for (uint64_t frame = 0; frame < 6; ++frame) {
        md_gpu_sync_wait(frame_done[frame % 2]);
        md_gpu_temp_t scope = md_gpu_temp_begin(f.gfx);
        md_gpu_render_desc_t rd = {0};
        rd.color_count = 1;
        rd.color[0].texture = target;
        rd.color[0].load    = MD_GPU_LOAD_CLEAR;
        rd.color[0].clear.f32[2] = 1.0f;      /* blue */
        rd.color[0].clear.f32[3] = 1.0f;
        rd.label = "ui";
        ASSERT_TRUE(md_gpu_render_begin(f.gfx, &rd));
        imgui_render(&ui, f.gfx, &dd);
        ASSERT_TRUE(md_gpu_render_end(f.gfx));
        md_gpu_temp_end(f.gfx, scope);
        frame_done[frame % 2] = md_gpu_stream_record(f.gfx);
        md_gpu_device_poll(f.dev);
    }

    static uint32_t px[32 * 32];
    ASSERT_TRUE(r_read(&f, target, 0, 0, px, sizeof(px)));
    EXPECT_EQ(rgba8(255, 255, 255, 255), px[4 * 32 + 4]);          /* white, inside its clip   */
    EXPECT_EQ(rgba8(0, 0, 255, 255),     px[24 * 32 + 4]);         /* clipped away: clear      */
    EXPECT_TRUE(near8(px[24 * 32 + 24], rgba8(128, 0, 127, 255), 1)); /* red at 50% over blue */

    md_gpu_pipeline_destroy(ui.pipeline);
    md_gpu_texture_destroy(target);
    md_gpu_texture_destroy(font);
    r_close(&f);
}

/* =========================================================================
   Presentation, through VK_EXT_headless_surface where the backend has it
   ========================================================================= */

static md_gpu_surface_t r_headless(rfix_t* f, uint32_t w, uint32_t h, md_gpu_tex_usage_t usage) {
    md_gpu_surface_desc_t sd = {0};
    sd.system = MD_GPU_WINDOW_HEADLESS;
    sd.width  = w;
    sd.height = h;
    sd.usage  = usage;
    sd.label  = "headless";
    return md_gpu_surface_create(f->dev, &sd);
}

UTEST(gpu_render, headless_surface_frames_resize_and_minimise) {
    rfix_t f;
    if (!r_open(&f)) UTEST_SKIP(r_skip_reason());
    if (!f.present) { r_close(&f); UTEST_SKIP("the device cannot present"); }
    md_gpu_surface_t sf = r_headless(&f, 32, 24, 0);
    if (!sf) { const char* why = r_skip_reason(); r_close(&f); UTEST_SKIP(why); }

    md_gpu_sync_t frame_done[2] = {md_gpu_sync_none(), md_gpu_sync_none()};
    md_gpu_mem_t rb = md_gpu_malloc(f.gfx, MD_GPU_MEM_HOST_READ, 48 * 40 * 4);
    ASSERT_TRUE(rb.cpu != NULL);
    for (uint64_t frame = 0; frame < 8; ++frame) {
        if (frame == 4) md_gpu_surface_resize(sf, 48, 40);
        md_gpu_sync_wait(frame_done[frame % 2]);
        md_gpu_texture_t back = md_gpu_surface_acquire(f.gfx, sf);
        ASSERT_TRUE(back != NULL);
        const md_gpu_texture_desc_t* d = md_gpu_texture_desc(back);
        EXPECT_EQ(frame < 4 ? 32u : 48u, d->width);
        EXPECT_EQ(frame < 4 ? 24u : 40u, d->height);
        EXPECT_EQ(MD_GPU_FORMAT_BGRA8_UNORM, d->format);
        EXPECT_FALSE(md_gpu_surface_acquire(f.gfx, sf) != NULL);      /* one at a time */

        md_gpu_render_desc_t rd = {0};
        rd.color_count = 1;
        rd.color[0].texture = back;
        rd.color[0].load    = MD_GPU_LOAD_CLEAR;
        rd.color[0].clear.f32[0] = (float)(frame % 4) / 3.0f;
        rd.color[0].clear.f32[3] = 1.0f;
        ASSERT_TRUE(md_gpu_render_begin(f.gfx, &rd));
        EXPECT_FALSE(md_gpu_surface_present(f.gfx, sf));             /* not inside a pass */
        ASSERT_TRUE(md_gpu_render_end(f.gfx));

        md_gpu_tex_region_t px = {0};
        px.offset[0] = d->width - 1; px.offset[1] = d->height - 1; px.extent[0] = px.extent[1] = px.extent[2] = 1;
        EXPECT_TRUE(md_gpu_copy_from_texture(f.gfx, rb.gpu, back, &px));
        ASSERT_TRUE(md_gpu_surface_present(f.gfx, sf));
        EXPECT_FALSE(md_gpu_surface_present(f.gfx, sf));             /* nothing acquired */
        frame_done[frame % 2] = md_gpu_stream_record(f.gfx);
        md_gpu_sync_wait(frame_done[frame % 2]);
        const uint8_t* bgra = (const uint8_t*)rb.cpu;
        const uint8_t r = (uint8_t)((float)(frame % 4) / 3.0f * 255.0f + 0.5f);
        EXPECT_NEAR((int)r, (int)bgra[2], 1);                        /* BGRA: red is byte 2 */
        EXPECT_EQ(255, (int)bgra[3]);
        md_gpu_device_poll(f.dev);
    }

    /* A zero-sized drawable (minimised window) skips the frame. */
    md_gpu_surface_resize(sf, 0, 0);
    EXPECT_TRUE(md_gpu_surface_acquire(f.gfx, sf) == NULL);
    md_gpu_surface_resize(sf, 16, 16);
    md_gpu_texture_t back = md_gpu_surface_acquire(f.gfx, sf);
    ASSERT_TRUE(back != NULL);
    EXPECT_EQ(16u, md_gpu_texture_desc(back)->width);
    ASSERT_TRUE(r_begin(&f, back, MD_GPU_LOAD_CLEAR, 0, 0, 0, 1, NULL, 0));
    ASSERT_TRUE(md_gpu_render_end(f.gfx));
    ASSERT_TRUE(md_gpu_surface_present(f.gfx, sf));

    md_gpu_free(f.gfx, rb.gpu);
    md_gpu_surface_destroy(sf);
    md_gpu_device_poll(f.dev);
    r_close(&f);
}

UTEST(gpu_render, surface_usage_and_format_requests) {
    rfix_t f;
    if (!r_open(&f)) UTEST_SKIP(r_skip_reason());
    if (!f.present) { r_close(&f); UTEST_SKIP("the device cannot present"); }

    /* A SAMPLED surface gets a sampled handle; a surface texture is not the
       caller's to destroy. */
    md_gpu_surface_t sf = r_headless(&f, 8, 8, MD_GPU_TEX_SAMPLED);
    if (sf) {
        md_gpu_texture_t back = md_gpu_surface_acquire(f.gfx, sf);
        ASSERT_TRUE(back != NULL);
        EXPECT_NE(0u, (unsigned)md_gpu_texture_sampled(back).handle);
        md_gpu_texture_destroy(back);            /* rejected */
        ASSERT_TRUE(r_begin(&f, back, MD_GPU_LOAD_CLEAR, 0, 0, 0, 1, NULL, 0));
        ASSERT_TRUE(md_gpu_render_end(f.gfx));
        ASSERT_TRUE(md_gpu_surface_present(f.gfx, sf));
        md_gpu_surface_destroy(sf);
    }
    /* A format the surface does not offer fails by name, with no substitute. */
    md_gpu_surface_desc_t sd = {0};
    sd.system = MD_GPU_WINDOW_HEADLESS;
    sd.width = sd.height = 8;
    sd.format = MD_GPU_FORMAT_R32_UINT;
    EXPECT_TRUE(md_gpu_surface_create(f.dev, &sd) == NULL);
    /* Acquire needs a GRAPHICS stream. */
    md_gpu_surface_t s2 = r_headless(&f, 8, 8, 0);
    if (s2) {
        EXPECT_TRUE(md_gpu_surface_acquire(md_gpu_stream_default(f.dev, MD_GPU_STREAM_COMPUTE), s2) == NULL);
        md_gpu_surface_destroy(s2);          /* never acquired: still fine */
    }
    r_close(&f);
}

#endif /* MD_ENABLE_GPU */
