#include "utest.h"

#include <md_render.h>
#include <md_render_mol.h>

#include <md_system.h>
#include <md_system_state.h>
#include <md_pdb.h>
#include <md_infer.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if MD_ENABLE_GPU

/* md_render against analytic answers: every covered pixel, depth and pick is
   compared with a ray intersection computed here on the CPU, for the same
   pixel-centre rays (sample 0 is unjittered). Pixels whose ray grazes a
   silhouette within a small margin are not judged. */

typedef struct {
    md_gpu_device_t  dev;
    md_gpu_stream_t  gfx;
    md_render_t      r;
    md_render_view_t view;
} fix_t;

static bool f_open(fix_t* f) {
    memset(f, 0, sizeof(*f));
    md_gpu_device_desc_t dd = {0};
    dd.enable_validation = true;
    dd.label = "md_render unittest";
    f->dev = md_gpu_device_create(&dd);
    if (!f->dev) return false;
    md_gpu_device_info_t info;
    md_gpu_device_info(f->dev, &info);
    if (!info.supports_graphics) {
        md_gpu_device_destroy(f->dev);
        f->dev = NULL;
        return false;
    }
    f->gfx  = md_gpu_stream_default(f->dev, MD_GPU_STREAM_GRAPHICS);
    f->r    = md_render_create(f->dev, NULL);
    f->view = f->r ? md_render_view_create(f->r, "test view") : NULL;
    return f->gfx && f->r && f->view;
}

static const char* f_skip_reason(void) {
    const char* err = md_gpu_last_error();
    return (err && err[0]) ? err : "No GPU device with graphics support";
}

static void f_close(fix_t* f) {
    if (f->view) md_render_view_destroy(f->gfx, f->view);
    if (f->r) md_render_destroy(f->r);
    md_gpu_stream_sync(f->gfx);
    md_gpu_device_destroy(f->dev);
}

static md_gpu_addr_t f_buf(fix_t* f, const void* data, size_t size) {
    md_gpu_addr_t a = md_gpu_malloc(f->gfx, MD_GPU_MEM_DEVICE, size).gpu;
    if (a && data) md_gpu_upload(f->gfx, a, data, size);
    return a;
}

static bool f_read(fix_t* f, md_gpu_texture_t t, void* dst, size_t size) {
    md_gpu_tex_region_t r = {0};
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

static bool f_pick(fix_t* f, uint32_t x, uint32_t y, md_render_pick_t* out) {
    if (!md_render_view_pick(f->gfx, f->view, x, y)) return false;
    md_gpu_stream_sync(f->gfx);
    return md_render_view_pick_result(f->view, out);
}

/* ---- Cameras ------------------------------------------------------------------ */

typedef struct {
    bool  ortho;
    float tan_half;     /* perspective: tan(fovy / 2) */
    float half_h;       /* ortho: half height of the view */
    float aspect;
    float n, far_;
} cam_t;

/* Reverse Z, view space looking down -Z; the view matrix is identity. */
static md_render_camera_t make_camera(const cam_t* c) {
    md_render_camera_t m;
    memset(&m, 0, sizeof(m));
    m.view[0] = m.view[5] = m.view[10] = m.view[15] = 1.0f;
    if (c->ortho) {
        m.proj[0]  = 1.0f / (c->half_h * c->aspect);
        m.proj[5]  = 1.0f / c->half_h;
        m.proj[10] = 1.0f / (c->far_ - c->n);
        m.proj[14] = c->far_ / (c->far_ - c->n);
        m.proj[15] = 1.0f;
    } else {                                    /* infinite far plane */
        m.proj[0]  = 1.0f / (c->tan_half * c->aspect);
        m.proj[5]  = 1.0f / c->tan_half;
        m.proj[11] = -1.0f;
        m.proj[14] = c->n;
    }
    return m;
}

typedef struct { double o[3], d[3]; } ray_t;

static ray_t pixel_ray(const cam_t* c, uint32_t w, uint32_t h, double px, double py) {
    double nx = px / w * 2.0 - 1.0, ny = 1.0 - py / h * 2.0;
    ray_t r;
    if (c->ortho) {
        r.o[0] = nx * c->half_h * c->aspect; r.o[1] = ny * c->half_h; r.o[2] = 0;
        r.d[0] = 0; r.d[1] = 0; r.d[2] = -1;
    } else {
        double dx = nx * c->tan_half * c->aspect, dy = ny * c->tan_half, l = sqrt(dx * dx + dy * dy + 1);
        r.o[0] = r.o[1] = r.o[2] = 0;
        r.d[0] = dx / l; r.d[1] = dy / l; r.d[2] = -1 / l;
    }
    return r;
}

static double dot3(const double* a, const double* b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

/* t of the hit and the ray's distance to the centre */
static double sphere_t(ray_t r, const double c[3], double rad, double* dist) {
    double oc[3] = { r.o[0] - c[0], r.o[1] - c[1], r.o[2] - c[2] };
    double b = dot3(oc, r.d), cc = dot3(oc, oc);
    double d2 = cc - b * b;
    *dist = sqrt(d2 > 0 ? d2 : 0);
    double h = rad * rad - d2;
    return h < 0 ? -1 : -b - sqrt(h);
}

/* Distance from the ray (line) to the segment, and t of the capsule hit (-1 none) */
static double capsule_t(ray_t r, const double a[3], const double b[3], double rad, double* dist) {
    /* closest approach by sampling the segment finely, then the exact hit by bisection on the ray */
    double best = 1e30;
    for (int i = 0; i <= 2000; ++i) {
        double s = i / 2000.0, p[3] = { a[0] + (b[0] - a[0]) * s, a[1] + (b[1] - a[1]) * s, a[2] + (b[2] - a[2]) * s };
        double sd;
        sphere_t(r, p, 0, &sd);
        if (sd < best) best = sd;
    }
    *dist = best;
    if (best >= rad) return -1;
    double tmin = 1e30;
    for (int i = 0; i <= 400; ++i) {
        double s = i / 400.0, p[3] = { a[0] + (b[0] - a[0]) * s, a[1] + (b[1] - a[1]) * s, a[2] + (b[2] - a[2]) * s };
        double sd, t = sphere_t(r, p, rad, &sd);
        if (t >= 0 && t < tmin) tmin = t;
    }
    /* the union of spheres is within (rad)(1 - cos) of the capsule: good to ~1e-5 at this sampling */
    return tmin < 1e29 ? tmin : -1;
}

static float proj_depth(const md_render_camera_t* m, const double p[3]) {
    double z = m->proj[10] * p[2] + m->proj[14];
    double w = m->proj[11] * p[2] + m->proj[15];
    return (float)(z / w);
}

/* ---- Scenes -------------------------------------------------------------------- */

/* sRGB colours go through linear and back in the resolve: within one step */
static bool near_rgba(uint32_t a, uint32_t b) {
    for (int c = 0; c < 4; ++c) {
        int x = (int)((a >> (8 * c)) & 0xFF), y = (int)((b >> (8 * c)) & 0xFF);
        if (abs(x - y) > 1) return false;
    }
    return true;
}

static md_render_frame_t frame_of(uint32_t w, uint32_t h, md_render_camera_t cam, const md_render_batch_t* b, uint32_t n) {
    md_render_frame_t f;
    memset(&f, 0, sizeof(f));
    f.width = w;
    f.height = h;
    f.camera = cam;
    f.batches = b;
    f.batch_count = n;
    f.jitter = MD_RENDER_JITTER_NONE;
    f.samples = 1;
    return f;
}

static int check_spheres(fix_t* fx, const cam_t* c, md_render_camera_t cam, uint32_t W, uint32_t H,
                         const double (*ctr)[3], const double* rad, uint32_t n, int* judged) {
    md_render_targets_t t = md_render_view_targets(fx->view);
    uint32_t* id = malloc(W * H * 8);
    float*    dp = malloc(W * H * 4);
    int bad = 0;
    *judged = 0;
    if (!f_read(fx, t.id, id, W * H * 8) || !f_read(fx, t.depth, dp, W * H * 4)) { free(id); free(dp); return -1; }
    for (uint32_t y = 0; y < H; ++y) {
        for (uint32_t x = 0; x < W; ++x) {
            ray_t r = pixel_ray(c, W, H, x + 0.5, y + 0.5);
            double best = 1e30; int who = -1; bool unsure = false;
            for (uint32_t i = 0; i < n; ++i) {
                double dist, tt = sphere_t(r, ctr[i], rad[i], &dist);
                if (fabs(dist - rad[i]) < 2e-3 * rad[i]) unsure = true;
                if (tt >= 0 && tt < best) { best = tt; who = (int)i; }
            }
            if (unsure) continue;
            (*judged)++;
            uint32_t ix = id[2 * (y * W + x)], iy = id[2 * (y * W + x) + 1];
            if (who < 0) {
                if (iy != 0) bad++;
                continue;
            }
            if ((iy >> 20) == 0 || ix != (uint32_t)who) { bad++; continue; }
            double p[3] = { r.o[0] + r.d[0] * best, r.o[1] + r.d[1] * best, r.o[2] + r.d[2] * best };
            float ed = proj_depth(&cam, p);
            if (fabsf(ed - dp[y * W + x]) > 1e-4f * fmaxf(fabsf(ed), 1e-3f) + 1e-6f) bad++;
        }
    }
    free(id); free(dp);
    return bad;
}

UTEST(render, sphere_coverage_depth_pick_perspective) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 96, H = 64;
    cam_t c = { .tan_half = 0.4f, .aspect = (float)W / H, .n = 0.1f };
    md_render_camera_t cam = make_camera(&c);

    const double ctr[3][3] = { { 0, 0, -10 }, { 2.2, 0.8, -12 }, { -1.5, -0.9, -8 } };
    const double rad[3]    = { 1.5, 2.0, 0.7 };
    float pos[9], radius[3];
    for (int i = 0; i < 3; ++i) { for (int k = 0; k < 3; ++k) pos[3 * i + k] = (float)ctr[i][k]; radius[i] = (float)rad[i]; }

    md_render_batch_t b = {0};
    b.key = 7;
    b.kind = MD_RENDER_SPHERES;
    b.count = 3;
    b.points.count = 3;
    b.points.position = f_buf(&fx, pos, sizeof(pos));
    b.points.radius = f_buf(&fx, radius, sizeof(radius));
    b.style.spheres.radius_scale = 1.0f;
    b.color.mode = MD_RENDER_COLOR_UNIFORM;
    b.color.uniform = 0xFF00FF00u;

    md_render_frame_t f = frame_of(W, H, cam, &b, 1);
    md_render_status_t st;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    EXPECT_EQ(MD_RENDER_WORK_RASTER | MD_RENDER_WORK_RESOLVE | MD_RENDER_WORK_MASK, st.work);

    int judged = 0;
    int bad = check_spheres(&fx, &c, cam, W, H, ctr, rad, 3, &judged);
    EXPECT_EQ(0, bad);
    EXPECT_GT(judged, (int)(W * H * 9 / 10));

    /* Pick the projected centre of sphere 2 */
    uint32_t px = (uint32_t)((ctr[2][0] / -ctr[2][2] / (c.tan_half * c.aspect) * 0.5 + 0.5) * W);
    uint32_t py = (uint32_t)((0.5 - ctr[2][1] / -ctr[2][2] / c.tan_half * 0.5) * H);
    md_render_pick_t pk;
    ASSERT_TRUE(f_pick(&fx, px, py, &pk));
    EXPECT_EQ(7u, (unsigned)pk.key);
    EXPECT_EQ(MD_RENDER_SPHERES, pk.kind);
    EXPECT_EQ(2u, pk.element);
    EXPECT_EQ(2u, pk.attr);
    EXPECT_EQ(2u, pk.id);
    ray_t r = pixel_ray(&c, W, H, px + 0.5, py + 0.5);
    double dist, tt = sphere_t(r, ctr[2], rad[2], &dist);
    for (int k = 0; k < 3; ++k) EXPECT_NEAR(r.o[k] + r.d[k] * tt, pk.position[k], 1e-3);

    /* Background */
    ASSERT_TRUE(f_pick(&fx, 0, 0, &pk));
    EXPECT_EQ(0u, (unsigned)pk.key);

    /* Colour: uniform green, covered */
    uint32_t* col = malloc(W * H * 4);
    ASSERT_TRUE(f_read(&fx, md_render_view_targets(fx.view).color, col, W * H * 4));
    EXPECT_TRUE(near_rgba(0xFF00FF00u, col[py * W + px]));
    EXPECT_EQ(0u, col[0]);
    free(col);
    f_close(&fx);
}

UTEST(render, sphere_coverage_orthographic) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 64, H = 64;
    cam_t c = { .ortho = true, .half_h = 3.0f, .aspect = 1.0f, .n = 0.0f, .far_ = 50.0f };
    md_render_camera_t cam = make_camera(&c);
    const double ctr[2][3] = { { 0.5, 0.2, -10 }, { -1.0, -1.0, -20 } };
    const double rad[2]    = { 1.2, 2.5 };
    float pos[6], radius[2];
    for (int i = 0; i < 2; ++i) { for (int k = 0; k < 3; ++k) pos[3 * i + k] = (float)ctr[i][k]; radius[i] = (float)rad[i]; }
    md_render_batch_t b = {0};
    b.key = 1; b.kind = MD_RENDER_SPHERES; b.count = 2; b.points.count = 2;
    b.points.position = f_buf(&fx, pos, sizeof(pos));
    b.points.radius   = f_buf(&fx, radius, sizeof(radius));
    b.style.spheres.radius_scale = 1.0f;
    md_render_frame_t f = frame_of(W, H, cam, &b, 1);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    int judged = 0;
    EXPECT_EQ(0, check_spheres(&fx, &c, cam, W, H, ctr, rad, 2, &judged));
    EXPECT_GT(judged, (int)(W * H * 9 / 10));
    f_close(&fx);
}

UTEST(render, capsule_coverage_depth_pick) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 80, H = 60;
    cam_t c = { .tan_half = 0.45f, .aspect = (float)W / H, .n = 0.1f };
    md_render_camera_t cam = make_camera(&c);
    /* Two capsules: one oblique, receding in depth; one nearly along the view axis */
    const double P[4][3] = { { -3, -1, -9 }, { 2.5, 1.5, -14 }, { 1.0, -1.5, -7 }, { 1.3, -1.2, -11 } };
    const double R = 0.6;
    float pos[12];
    for (int i = 0; i < 4; ++i) for (int k = 0; k < 3; ++k) pos[3 * i + k] = (float)P[i][k];
    uint32_t colors[4] = { 0xFF0000FFu, 0xFFFF0000u, 0xFF00FF00u, 0xFF00FFFFu };
    md_render_batch_t b = {0};
    b.key = 3; b.kind = MD_RENDER_CAPSULES; b.count = 2; b.points.count = 4;
    b.points.position = f_buf(&fx, pos, sizeof(pos));
    b.points.color    = f_buf(&fx, colors, sizeof(colors));
    b.style.capsules.radius = (float)R;
    b.color.sharpness = 1.0f;
    md_render_frame_t f = frame_of(W, H, cam, &b, 1);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));

    md_render_targets_t t = md_render_view_targets(fx.view);
    uint32_t* id = malloc(W * H * 8);
    float*    dp = malloc(W * H * 4);
    ASSERT_TRUE(f_read(&fx, t.id, id, W * H * 8));
    ASSERT_TRUE(f_read(&fx, t.depth, dp, W * H * 4));
    int bad = 0, judged = 0, covered = 0;
    for (uint32_t y = 0; y < H; ++y) for (uint32_t x = 0; x < W; ++x) {
        ray_t r = pixel_ray(&c, W, H, x + 0.5, y + 0.5);
        double best = 1e30; int who = -1; bool unsure = false;
        for (int e = 0; e < 2; ++e) {
            double dist, tt = capsule_t(r, P[2 * e], P[2 * e + 1], R, &dist);
            if (fabs(dist - R) < 3e-3 * R) unsure = true;
            if (tt >= 0 && tt < best) { best = tt; who = e; }
        }
        if (unsure) continue;
        judged++;
        uint32_t ix = id[2 * (y * W + x)], iy = id[2 * (y * W + x) + 1];
        if (who < 0) { if (iy) bad++; continue; }
        covered++;
        if (!iy || ix != (uint32_t)who) { bad++; continue; }
        double p[3] = { r.d[0] * best, r.d[1] * best, r.d[2] * best };
        float ed = proj_depth(&cam, p);
        if (fabsf(ed - dp[y * W + x]) > 2e-4f * ed) bad++;
    }
    free(id); free(dp);
    EXPECT_EQ(0, bad);
    EXPECT_GT(covered, 200);
    EXPECT_GT(judged, (int)(W * H * 9 / 10));

    /* Pick near each end of capsule 0: the nearer point's attribute, its colour (sharp) */
    for (int end = 0; end < 2; ++end) {
        double s = end ? 0.8 : 0.2, q[3];
        for (int k = 0; k < 3; ++k) q[k] = P[0][k] + (P[1][k] - P[0][k]) * s;
        uint32_t px = (uint32_t)((q[0] / -q[2] / (c.tan_half * c.aspect) * 0.5 + 0.5) * W);
        uint32_t py = (uint32_t)((0.5 - q[1] / -q[2] / c.tan_half * 0.5) * H);
        md_render_pick_t pk;
        ASSERT_TRUE(f_pick(&fx, px, py, &pk));
        EXPECT_EQ(3u, (unsigned)pk.key);
        EXPECT_EQ(0u, pk.element);
        EXPECT_EQ((uint32_t)end, pk.attr);
        uint32_t* col = malloc(W * H * 4);
        ASSERT_TRUE(f_read(&fx, t.color, col, W * H * 4));
        EXPECT_TRUE(near_rgba(colors[end], col[py * W + px]));
        free(col);
    }
    f_close(&fx);
}

UTEST(render, index_visibility_instances) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 64, H = 64;
    cam_t c = { .ortho = true, .half_h = 4.0f, .aspect = 1.0f, .n = 0.0f, .far_ = 100.0f };
    md_render_camera_t cam = make_camera(&c);
    /* 4 points on a row; draw points 3 and 1 through an index; hide point 1 */
    float pos[12] = { -3, 0, -10,  -1, 0, -10,  1, 0, -10,  3, 0, -10 };
    uint32_t index[2] = { 3, 1 };
    uint32_t visible = ~(1u << 1);
    md_render_batch_t b = {0};
    b.key = 11; b.kind = MD_RENDER_SPHERES; b.count = 2; b.points.count = 4;
    b.points.position   = f_buf(&fx, pos, sizeof(pos));
    b.points.visible[1] = f_buf(&fx, &visible, 4);
    b.index = f_buf(&fx, index, sizeof(index));
    b.style.spheres.radius_add = 0.8f;
    /* Two instances: as is, and moved down by 2 */
    float xf[32] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1,   1,0,0,0, 0,1,0,0, 0,0,1,0, 0,-2,0,1 };
    b.transform = xf;
    b.transform_count = 2;
    md_render_frame_t f = frame_of(W, H, cam, &b, 1);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));

    /* world (x, y) -> pixel */
    #define PX(wx) ((uint32_t)(((wx) / 4.0 * 0.5 + 0.5) * W))
    #define PY(wy) ((uint32_t)((0.5 - (wy) / 4.0 * 0.5) * H))
    md_render_pick_t pk;
    ASSERT_TRUE(f_pick(&fx, PX(3.0), PY(0.0), &pk));
    EXPECT_EQ(11u, (unsigned)pk.key);
    EXPECT_EQ(0u, pk.element);          /* element 0 is point 3 */
    EXPECT_EQ(3u, pk.attr);
    EXPECT_EQ(0u, pk.instance);
    ASSERT_TRUE(f_pick(&fx, PX(3.0), PY(-2.0), &pk));
    EXPECT_EQ(11u, (unsigned)pk.key);
    EXPECT_EQ(1u, pk.instance);
    EXPECT_NEAR(4.0 * (1.0 - 2.0 * (PY(-2.0) + 0.5) / H), pk.position[1], 1e-4);   /* the pixel centre's ray */
    ASSERT_TRUE(f_pick(&fx, PX(-1.0), PY(0.0), &pk));   /* hidden */
    EXPECT_EQ(0u, (unsigned)pk.key);
    ASSERT_TRUE(f_pick(&fx, PX(1.0), PY(0.0), &pk));    /* not indexed */
    EXPECT_EQ(0u, (unsigned)pk.key);
    #undef PX
    #undef PY
    f_close(&fx);
}

UTEST(render, lazy_tiers) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 32, H = 32;
    cam_t c = { .tan_half = 0.5f, .aspect = 1.0f, .n = 0.1f };
    float pos[3] = { 0, 0, -5 };
    uint32_t sel = 1;
    md_render_batch_t b = {0};
    b.key = 1; b.kind = MD_RENDER_SPHERES; b.count = 1; b.points.count = 1;
    b.points.position = f_buf(&fx, pos, sizeof(pos));
    b.points.state[0] = f_buf(&fx, &sel, 4);
    b.style.spheres.radius_add = 1.0f;
    b.color.mode = MD_RENDER_COLOR_UNIFORM;
    b.color.uniform = 0xFFFFFFFFu;
    md_render_frame_t f = frame_of(W, H, make_camera(&c), &b, 1);
    f.jitter = MD_RENDER_JITTER_STATIC;
    f.samples = 4;

    md_render_status_t st;
    const md_render_work_t ALL = MD_RENDER_WORK_RASTER | MD_RENDER_WORK_RESOLVE | MD_RENDER_WORK_MASK;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    EXPECT_EQ(ALL, st.work);
    EXPECT_EQ(0u, st.sample);
    EXPECT_EQ(0.0f, st.jitter[0]);
    /* Still: the jitter sequence runs to `samples`, then nothing is drawn */
    for (uint32_t k = 1; k < 4; ++k) {
        ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
        EXPECT_EQ(ALL, st.work);
        EXPECT_EQ(k, st.sample);
        EXPECT_TRUE(st.jitter[0] != 0.0f || st.jitter[1] != 0.0f);
    }
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    EXPECT_EQ(0u, st.work);
    EXPECT_EQ(3u, st.sample);

    /* State: masks only, the sequence carries on */
    f.content[MD_RENDER_CONTENT_STATE]++;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    EXPECT_EQ(MD_RENDER_WORK_MASK, st.work);
    EXPECT_EQ(3u, st.sample);

    /* Appearance, through the frame (a colour field): resolve, sample 0 */
    b.color.uniform = 0xFF0000FFu;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    EXPECT_EQ(MD_RENDER_WORK_RESOLVE | MD_RENDER_WORK_MASK, st.work);
    EXPECT_EQ(0u, st.sample);
    /* ...and the next draw continues the sequence with a raster */
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    EXPECT_EQ(ALL, st.work);
    EXPECT_EQ(1u, st.sample);

    /* Appearance, through the counter */
    f.content[MD_RENDER_CONTENT_APPEARANCE]++;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    EXPECT_EQ(MD_RENDER_WORK_RESOLVE | MD_RENDER_WORK_MASK, st.work);

    /* Geometry, through the frame (a style) and through the counter */
    for (uint32_t k = 0; k < 8; ++k) ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    EXPECT_EQ(0u, st.work);
    b.style.spheres.radius_add = 1.1f;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    EXPECT_EQ(ALL, st.work);
    EXPECT_EQ(0u, st.sample);
    for (uint32_t k = 0; k < 8; ++k) ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    f.content[MD_RENDER_CONTENT_GEOMETRY]++;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    EXPECT_EQ(ALL, st.work);

    /* Camera */
    for (uint32_t k = 0; k < 8; ++k) ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    f.camera.view[14] = -0.5f;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    EXPECT_EQ(ALL, st.work);

    /* Invalid frames change nothing */
    md_render_batch_t bad = b;
    bad.key = 0;
    f.batches = &bad;
    EXPECT_FALSE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    f_close(&fx);
}

UTEST(render, color_modes) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 64, H = 16;
    cam_t c = { .ortho = true, .half_h = 1.0f, .aspect = 4.0f, .n = 0.0f, .far_ = 100.0f };
    float pos[12] = { -3, 0, -10,  -1, 0, -10,  1, 0, -10,  3, 0, -10 };
    uint32_t colors[4]  = { 0xFF0000FFu, 0xFF00FF00u, 0xFFFF0000u, 0xFF808080u };
    uint16_t pal_idx[4] = { 2, 0, 1, 7 };        /* 7: past the palette, takes entry 0 */
    uint32_t palette[3] = { 0xFF102030u, 0xFF405060u, 0xFF708090u };
    md_render_batch_t b[2] = {{0}};
    b[0].key = 1; b[0].kind = MD_RENDER_SPHERES; b[0].count = 2; b[0].points.count = 4;
    b[0].points.position = f_buf(&fx, pos, sizeof(pos));
    b[0].points.color    = f_buf(&fx, colors, sizeof(colors));
    b[0].style.spheres.radius_add = 0.5f;
    b[1] = b[0];
    b[1].key = 2;
    uint32_t idx2[2] = { 2, 3 };
    b[1].index = f_buf(&fx, idx2, sizeof(idx2));
    b[1].points.palette_index = f_buf(&fx, pal_idx, sizeof(pal_idx));
    b[1].color.mode = MD_RENDER_COLOR_PALETTE;
    b[1].color.palette = f_buf(&fx, palette, sizeof(palette));
    b[1].color.palette_count = 3;
    md_render_frame_t f = frame_of(W, H, make_camera(&c), b, 2);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    uint32_t* col = malloc(W * H * 4);
    ASSERT_TRUE(f_read(&fx, md_render_view_targets(fx.view).color, col, W * H * 4));
    #define AT(wx) col[(H / 2) * W + (uint32_t)(((wx) / 4.0 * 0.5 + 0.5) * W)]
    EXPECT_TRUE(near_rgba(colors[0], AT(-3.0)));
    EXPECT_TRUE(near_rgba(colors[1], AT(-1.0)));
    EXPECT_TRUE(near_rgba(palette[1], AT(1.0)));
    EXPECT_TRUE(near_rgba(palette[0], AT(3.0)));

    /* Tint fully to black, desaturate to grey: appearance only */
    b[0].color.tint[3] = 1.0f;     /* tint rgb 0 */
    b[1].color.desaturate = 1.0f;
    md_render_status_t st;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    EXPECT_EQ(MD_RENDER_WORK_RESOLVE | MD_RENDER_WORK_MASK, st.work);
    ASSERT_TRUE(f_read(&fx, md_render_view_targets(fx.view).color, col, W * H * 4));
    EXPECT_EQ(0xFF000000u, AT(-3.0));
    uint32_t g = AT(1.0);
    EXPECT_EQ(g & 0xFF, (g >> 8) & 0xFF);
    EXPECT_EQ(g & 0xFF, (g >> 16) & 0xFF);
    #undef AT
    free(col);
    f_close(&fx);
}

UTEST(render, motion_from_velocity) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 64, H = 64;
    cam_t c = { .tan_half = 0.5f, .aspect = 1.0f, .n = 0.1f };
    md_render_camera_t cam = make_camera(&c);
    float pos[3] = { 0, 0, -10 };
    float vel[3] = { 2.0f, -1.0f, 0.0f };      /* units per second */
    uint32_t pv[2];
    md_render_pack_velocity(vel, pv);
    md_gpu_addr_t P = f_buf(&fx, pos, sizeof(pos));
    md_render_batch_t b = {0};
    b.key = 1; b.kind = MD_RENDER_SPHERES; b.count = 1; b.points.count = 1;
    b.points.position = P;
    b.points.velocity = f_buf(&fx, pv, sizeof(pv));
    b.style.spheres.radius_add = 1.0f;
    md_render_frame_t f = frame_of(W, H, cam, &b, 1);
    f.time = 0.0;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));

    /* 0.25 s later the point has moved by v * dt; its motion is that, projected */
    const double dt = 0.25;
    float pos2[3] = { (float)(pos[0] + vel[0] * dt), (float)(pos[1] + vel[1] * dt), pos[2] };
    md_gpu_upload(fx.gfx, P, pos2, sizeof(pos2));
    f.content[MD_RENDER_CONTENT_GEOMETRY]++;
    f.time = dt;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));

    float* mo = malloc(W * H * 4);
    uint16_t* raw = (uint16_t*)mo;
    ASSERT_TRUE(f_read(&fx, md_render_view_targets(fx.view).motion, mo, W * H * 4));
    /* The sphere's front point under its new centre's pixel */
    double scale = 1.0 / (c.tan_half * 10.0) * 0.5 * W;   /* pixels per unit at depth 10 (front at 9: close) */
    uint32_t px = (uint32_t)(pos2[0] / 10.0 / c.tan_half * 0.5 * W + W / 2);
    uint32_t py = (uint32_t)(H / 2 - pos2[1] / 10.0 / c.tan_half * 0.5 * H);
    /* half floats */
    uint16_t hx = raw[2 * (py * W + px)], hy = raw[2 * (py * W + px) + 1];
    float mx, my;
    {
        uint32_t s, e, m, bits;
        #define H2F(h, out) s = (h >> 15) & 1; e = (h >> 10) & 31; m = h & 1023; \
            bits = e ? ((s << 31) | ((e + 112) << 23) | (m << 13)) : (s << 31); memcpy(&out, &bits, 4);
        H2F(hx, mx) H2F(hy, my)
        #undef H2F
    }
    /* The hit point in front of the centre is at depth ~9, so the projected motion is v dt / 9 rather than / 10 */
    double ex = vel[0] * dt / 9.0 / c.tan_half * 0.5 * W;
    double ey = -vel[1] * dt / 9.0 / c.tan_half * 0.5 * H;
    (void)scale;
    EXPECT_NEAR(ex, mx, 0.15 * fabs(ex));
    EXPECT_NEAR(ey, my, 0.15 * fabs(ey));
    /* Background: no motion */
    EXPECT_EQ(0, raw[0]);
    free(mo);
    f_close(&fx);
}

UTEST(render, masks_visible_and_hidden) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 64, H = 32;
    cam_t c = { .ortho = true, .half_h = 2.0f, .aspect = 2.0f, .n = 0.0f, .far_ = 100.0f };
    /* Point 0 selected, behind point 1 which overlaps its right half; point 2 highlighted */
    float pos[9] = { -1.0f, 0, -20,   0.0f, 0, -10,   2.5f, 0, -10 };
    uint32_t sel = 1u << 0, hl = 1u << 2;
    md_render_batch_t b = {0};
    b.key = 1; b.kind = MD_RENDER_SPHERES; b.count = 3; b.points.count = 3;
    b.points.position = f_buf(&fx, pos, sizeof(pos));
    b.points.state[0] = f_buf(&fx, &sel, 4);
    b.points.state[1] = f_buf(&fx, &hl, 4);
    b.style.spheres.radius_add = 0.9f;
    md_render_frame_t f = frame_of(W, H, make_camera(&c), &b, 1);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    uint32_t* m = malloc(W * H * 4);
    ASSERT_TRUE(f_read(&fx, md_render_view_targets(fx.view).mask, m, W * H * 4));
    #define AT(wx) m[(H / 2) * W + (uint32_t)(((wx) / 4.0 * 0.5 + 0.5) * W)]
    EXPECT_EQ(0x0000FF00u | 0x000000FFu, AT(-1.5));   /* selection, visible and there   */
    EXPECT_EQ(0x0000FF00u, AT(-0.3));                 /* selection hidden behind point 1 */
    EXPECT_EQ(0xFF000000u | 0x00FF0000u, AT(2.5));    /* highlight, visible             */
    EXPECT_EQ(0u, AT(0.6));                           /* point 1: nothing               */
    #undef AT

    /* Selection moves to point 1: STATE only */
    sel = 1u << 1;
    md_gpu_upload(fx.gfx, b.points.state[0], &sel, 4);
    f.content[MD_RENDER_CONTENT_STATE]++;
    md_render_status_t st;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    EXPECT_EQ(MD_RENDER_WORK_MASK, st.work);
    ASSERT_TRUE(f_read(&fx, md_render_view_targets(fx.view).mask, m, W * H * 4));
    #define AT(wx) m[(H / 2) * W + (uint32_t)(((wx) / 4.0 * 0.5 + 0.5) * W)]
    EXPECT_EQ(0u, AT(-1.5));
    EXPECT_EQ(0x0000FF00u | 0x000000FFu, AT(-0.3));
    #undef AT
    free(m);
    f_close(&fx);
}

UTEST(render, external_ids) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    md_render_batch_t b[2] = {{0}};
    float pos[3] = { 0, 0, -5 };
    b[0].key = 5; b[0].kind = MD_RENDER_SPHERES; b[0].count = 1; b[0].points.count = 1;
    b[0].points.position = f_buf(&fx, pos, sizeof(pos));
    b[0].style.spheres.radius_add = 0.1f;
    b[1].key = 99; b[1].kind = MD_RENDER_EXTERNAL; b[1].count = 1000;
    cam_t c = { .tan_half = 0.5f, .aspect = 1.0f, .n = 0.1f };
    md_render_frame_t f = frame_of(16, 16, make_camera(&c), b, 2);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    /* The caller's renderer writes id 1 + 3 * 1000 + 17 (instance 3, element 17) at pixel (1, 1) */
    uint32_t id[2];
    md_render_encode_id(1, 0, 3 * 1000 + 17, id);
    md_render_targets_t t = md_render_view_targets(fx.view);
    md_gpu_mem_t up = md_gpu_malloc(fx.gfx, MD_GPU_MEM_HOST_WRITE, 8);
    memcpy(up.cpu, id, 8);
    md_gpu_tex_region_t reg = {0};
    reg.offset[0] = 1; reg.offset[1] = 1;
    reg.extent[0] = 1; reg.extent[1] = 1; reg.extent[2] = 1;
    ASSERT_TRUE(md_gpu_copy_to_texture(fx.gfx, t.id, &reg, up.gpu));
    md_render_pick_t pk;
    ASSERT_TRUE(f_pick(&fx, 1, 1, &pk));
    EXPECT_EQ(99u, (unsigned)pk.key);
    EXPECT_EQ(MD_RENDER_EXTERNAL, pk.kind);
    EXPECT_EQ(3u, pk.instance);
    EXPECT_EQ(17u, pk.element);
    EXPECT_EQ(UINT32_MAX, pk.attr);
    md_gpu_free(fx.gfx, up.gpu);
    f_close(&fx);
}

/* ---- The molecular layer ------------------------------------------------------- */

static const char ethanol_pdb[] =
    "HETATM    1  C1  EOH A   1      -0.748   0.015   0.024  1.00  0.00           C\n"
    "HETATM    2  C2  EOH A   1       0.558  -0.743  -0.059  1.00  0.00           C\n"
    "HETATM    3  O   EOH A   1       1.578   0.215   0.208  1.00  0.00           O\n"
    "HETATM    4  H11 EOH A   1      -0.832   0.507   0.996  1.00  0.00           H\n"
    "HETATM    5  H12 EOH A   1      -0.829   0.770  -0.762  1.00  0.00           H\n"
    "HETATM    6  H13 EOH A   1      -1.587  -0.676  -0.081  1.00  0.00           H\n"
    "HETATM    7  H21 EOH A   1       0.646  -1.239  -1.030  1.00  0.00           H\n"
    "HETATM    8  H22 EOH A   1       0.630  -1.501   0.726  1.00  0.00           H\n"
    "HETATM    9  HO  EOH A   1       2.414  -0.237   0.135  1.00  0.00           H\n"
    "END\n";

UTEST(render, mol_licorice_emit_and_pick) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, md_str_from_cstr(ethanol_pdb), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE, NULL));
    ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
    ASSERT_EQ(8u, (uint32_t)(sys.bond.count + sys.bond.noncovalent_count));

    uint64_t content[MD_RENDER_CONTENT_COUNT] = {0};

    /* Topology: plan, the caller allocates, write */
    md_temp_scope_t temp = md_temp_begin();
    md_render_mol_plan_t plan;
    ASSERT_TRUE(md_render_mol_plan(&plan, &sys, temp));
    md_gpu_mem_t mem = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, plan.size);
    md_render_mol_topo_t topo;
    ASSERT_TRUE(md_render_mol_write(fx.gfx, &topo, &plan, mem.gpu, content));
    md_temp_end(temp);
    EXPECT_EQ(sys.topology_version, topo.topology_version);
    EXPECT_EQ(8u, topo.bond_count);
    EXPECT_EQ(8u, topo.bond_class_offset[1]);   /* all covalent */

    /* The caller's atom arrays */
    md_gpu_addr_t P = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, 9 * 12).gpu;
    ASSERT_TRUE(md_render_mol_upload_state(fx.gfx, P, &state, content));
    uint32_t colors[9];
    for (int i = 0; i < 9; ++i) colors[i] = 0xFF000000u | (uint32_t)(i * 20);
    md_render_mol_atoms_t atoms = {0};
    atoms.position = P;

    md_render_mol_rep_t rep = {0};
    rep.key  = 42;
    rep.type = MD_RENDER_MOL_REP_LICORICE;
    rep.params.licorice.radius = 0.15f;
    rep.atom_color = f_buf(&fx, colors, sizeof(colors));
    rep.color.mode = MD_RENDER_COLOR_POINTS;

    md_render_batch_t batches[8];
    uint32_t nb = md_render_mol_emit(&topo, &atoms, &rep, batches, 8);
    ASSERT_EQ(2u, nb);                          /* atoms, covalent bonds */
    EXPECT_EQ(MD_RENDER_CAPSULES, batches[0].kind);
    EXPECT_EQ(MD_RENDER_SPHERES, batches[1].kind);
    EXPECT_EQ(8u, batches[0].count);
    /* emit is pure: asking with no room says how many it needs */
    EXPECT_EQ(2u, md_render_mol_emit(&topo, &atoms, &rep, NULL, 0));

    /* Looking down -Z at the molecule from 6 Å, orthographic */
    const uint32_t W = 128, H = 96;
    cam_t c = { .ortho = true, .half_h = 2.0f, .aspect = (float)W / H, .n = 0.0f, .far_ = 100.0f };
    md_render_camera_t cam = make_camera(&c);
    cam.view[14] = -6.0f;
    md_render_frame_t f = frame_of(W, H, cam, batches, nb);
    memcpy(f.content, content, sizeof(content));
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));

    /* The oxygen (atom 2) at its projected centre is an atom hit... */
    const float* o = &state.xyz[2].x;
    uint32_t px = (uint32_t)((o[0] / (c.half_h * c.aspect) * 0.5 + 0.5) * W);
    uint32_t py = (uint32_t)((0.5 - o[1] / c.half_h * 0.5) * H);
    md_render_pick_t pk;
    md_render_mol_hit_t hit;
    ASSERT_TRUE(f_pick(&fx, px, py, &pk));
    ASSERT_TRUE(md_render_mol_decode_pick(&pk, &hit));
    EXPECT_EQ(MD_RENDER_MOL_HIT_ATOM, hit.kind);
    EXPECT_EQ(42u, (unsigned)hit.rep_key);
    EXPECT_EQ(2u, hit.atom);

    /* ...and the middle of the C-O bond is a bond hit on that bond */
    const float* c2 = &state.xyz[1].x;
    float mid[2] = { 0.5f * (o[0] + c2[0]), 0.5f * (o[1] + c2[1]) };
    px = (uint32_t)((mid[0] / (c.half_h * c.aspect) * 0.5 + 0.5) * W);
    py = (uint32_t)((0.5 - mid[1] / c.half_h * 0.5) * H);
    ASSERT_TRUE(f_pick(&fx, px, py, &pk));
    ASSERT_TRUE(md_render_mol_decode_pick(&pk, &hit));
    EXPECT_EQ(MD_RENDER_MOL_HIT_BOND, hit.kind);
    ASSERT_LT(hit.bond, 8u);
    md_atom_pair_t pr = sys.bond.pairs[hit.bond];
    EXPECT_TRUE((pr.idx[0] == 1 && pr.idx[1] == 2) || (pr.idx[0] == 2 && pr.idx[1] == 1));

    if (getenv("MD_RENDER_DUMP")) {
        uint32_t* col = malloc(W * H * 4);
        if (f_read(&fx, md_render_view_targets(fx.view).color, col, W * H * 4)) {
            FILE* fp = fopen(getenv("MD_RENDER_DUMP"), "wb");
            if (fp) {
                fprintf(fp, "P6\n%u %u\n255\n", W, H);
                for (uint32_t i = 0; i < W * H; ++i) { uint8_t px3[3] = { (uint8_t)col[i], (uint8_t)(col[i] >> 8), (uint8_t)(col[i] >> 16) }; fwrite(px3, 1, 3, fp); }
                fclose(fp);
            }
        }
        free(col);
    }

    md_gpu_free(fx.gfx, P);
    md_gpu_free(fx.gfx, topo.base);
    md_system_free(&sys);
    f_close(&fx);
}

#endif
