#include "utest.h"

#include <md_render.h>
#include <md_render_mol.h>

#include <md_system.h>
#include <md_system_state.h>
#include <md_unitcell.h>
#include <md_pdb.h>
#include <md_gro.h>
#include <md_infer.h>
#include <md_backbone.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bench_rev.h"                    /* MD_BENCH_SOURCE_REV, generated at build time */

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
    static bool rev_printed = false;            /* once per run: a stale build shows as one */
    if (!rev_printed) { printf("    mdlib %s\n", MD_BENCH_SOURCE_REV); rev_printed = true; }
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

/* Radii from a table: float4 points whose w is a type (low 16 bits) with
   flags above it; the spheres as drawn are those of the table's radii, and a
   type past the table's end is radius 0 (not drawn) */
UTEST(render, sphere_radius_table) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 64, H = 64;
    cam_t c = { .ortho = true, .half_h = 3.0f, .aspect = 1.0f, .n = 0.0f, .far_ = 50.0f };
    md_render_camera_t cam = make_camera(&c);
    const double ctr[2][3] = { { 0.5, 0.2, -10 }, { -1.0, -1.0, -20 } };
    const double rad[2]    = { 1.2, 2.5 };
    const float table[3]   = { 0.0f, 2.5f, 1.2f };
    const uint32_t type[3] = { 2u | (0xABCDu << 16), 1u | (0x1234u << 16), 7u };   /* the third: past the table */
    float pts[12];
    for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < 3; ++k) pts[4 * i + k] = i < 2 ? (float)ctr[i][k] : (k == 2 ? -10.0f : 2.0f);
        memcpy(&pts[4 * i + 3], &type[i], 4);
    }
    md_gpu_addr_t P = f_buf(&fx, pts, sizeof(pts));
    md_render_batch_t b = {0};
    b.key = 1; b.kind = MD_RENDER_SPHERES; b.count = 3; b.points.count = 3;
    b.points.position = P;      b.points.stride = 16;
    b.points.radius   = P + 12; b.points.radius_stride = 16;
    b.points.radius_table = f_buf(&fx, table, sizeof(table));
    b.points.radius_table_count = 3;
    b.style.spheres.radius_scale = 1.0f;
    md_render_frame_t f = frame_of(W, H, cam, &b, 1);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    int judged = 0;
    EXPECT_EQ(0, check_spheres(&fx, &c, cam, W, H, ctr, rad, 2, &judged));
    EXPECT_GT(judged, (int)(W * H * 9 / 10));
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

/* ---- Splines ------------------------------------------------------------------- */

static int compare_views(fix_t* fx, md_render_view_t a, md_render_view_t b, uint32_t W, uint32_t H, int* covered);

/* Distance from p to the rectangle x in [x0, x1], y in [yc - a, yc + a], z 0 */
static double rect_dist(const double p[3], double x0, double x1, double yc, double a) {
    double dx = p[0] < x0 ? x0 - p[0] : p[0] > x1 ? p[0] - x1 : 0.0;
    double dy = fabs(p[1] - yc) - a;
    if (dy < 0) dy = 0;
    return sqrt(dx * dx + dy * dy + p[2] * p[2]);
}

/* Straight splines along local x: each part is a rectangle (the swept
   segment of its width) widened by the radius, exactly what the pieces of a
   straight spline make, caps included */
typedef struct bar_part_t { double x0, x1, yc, a; int element; } bar_part_t;
typedef struct bar_scene_t {
    double R[3][3];      /* local to view: view = R local + t */
    double t[3];
    double radius;
    const bar_part_t* part;
    int part_count;
} bar_scene_t;

/* The scene's distance at view-space p, and the nearest part and how far the next nearest is behind it */
static double bar_sdf(const bar_scene_t* s, const double pv[3], int* who, double* gap) {
    double q[3] = { pv[0] - s->t[0], pv[1] - s->t[1], pv[2] - s->t[2] }, p[3];
    for (int i = 0; i < 3; ++i) p[i] = s->R[0][i] * q[0] + s->R[1][i] * q[1] + s->R[2][i] * q[2];   /* R^T q */
    double best = 1e30, second = 1e30;
    *who = -1;
    for (int k = 0; k < s->part_count; ++k) {
        double d = rect_dist(p, s->part[k].x0, s->part[k].x1, s->part[k].yc, s->part[k].a) - s->radius;
        if (d < best) { second = best; best = d; *who = s->part[k].element; }
        else if (d < second) second = d;
    }
    *gap = second - best;
    return best;
}

/* Every pixel against the scene, sphere traced: coverage, the part (where no
   other is as near), depth. Rays that graze it are not judged. */
static int check_bars(fix_t* fx, const cam_t* c, md_render_camera_t cam, uint32_t W, uint32_t H, const bar_scene_t* s, int* judged) {
    md_render_targets_t t = md_render_view_targets(fx->view);
    uint32_t* id = malloc(W * H * 8);
    float*    dp = malloc(W * H * 4);
    int bad = 0;
    *judged = 0;
    if (!f_read(fx, t.id, id, W * H * 8) || !f_read(fx, t.depth, dp, W * H * 4)) { free(id); free(dp); return -1; }
    for (uint32_t y = 0; y < H; ++y) {
        for (uint32_t x = 0; x < W; ++x) {
            ray_t r = pixel_ray(c, W, H, x + 0.5, y + 0.5);
            double tt = 0.0, dmin = 1e30, gap = 0.0;
            int who = -1;
            bool hit = false;
            for (int it = 0; it < 4000 && tt < 200.0; ++it) {
                double p[3] = { r.o[0] + r.d[0] * tt, r.o[1] + r.d[1] * tt, r.o[2] + r.d[2] * tt };
                double d = bar_sdf(s, p, &who, &gap);
                if (d < 1e-7) { hit = true; break; }
                if (d < dmin) dmin = d;
                tt += d;
            }
            if (!hit && dmin < 2e-3) continue;              /* grazes */
            if (hit && dmin < 2e-3 && tt > 0) {
                /* near the silhouette on the way in: judge only if the ray came in steeply */
                continue;
            }
            (*judged)++;
            uint32_t ix = id[2 * (y * W + x)], iy = id[2 * (y * W + x) + 1];
            if (!hit) { if (iy != 0) bad++; continue; }
            if ((iy >> 20) == 0) { bad++; continue; }
            if (gap > 1e-3 && ix != (uint32_t)who) { bad++; if (bad < 6) printf("    pixel (%u, %u): part %u, expected %d\n", x, y, ix, who); continue; }
            double p[3] = { r.o[0] + r.d[0] * tt, r.o[1] + r.d[1] * tt, r.o[2] + r.d[2] * tt };
            float ed = proj_depth(&cam, p);
            if (fabsf(ed - dp[y * W + x]) > 1e-4f * fmaxf(fabsf(ed), 1e-3f) + 1e-6f) {
                bad++;
                if (bad < 6) printf("    pixel (%u, %u): depth %.7g, expected %.7g\n", x, y, dp[y * W + x], ed);
            }
        }
    }
    free(id); free(dp);
    return bad;
}

/* Straight splines are exact: a ribbon of seven control points whose parts
   are widths of their own (steps: caps on the wider side), one part hidden
   (caps on its neighbours), and a tube of four (width 0: a capsule). Turned,
   seen in perspective; coverage, the part each pixel shows and depth against
   the rectangles widened by the radius. Then hidden parts, a pick. */
UTEST(render, spline_straight_parts) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 128, H = 96;
    cam_t c = { .tan_half = 0.45f, .aspect = (float)W / H, .n = 0.1f };
    md_render_camera_t cam = make_camera(&c);
    enum { NA = 7, NB = 4, N = NA + NB };
    const double d = 1.2, radius = 0.22;
    const double wa[NA] = { 0.5, 0.5, 0.9, 0.9, 0.3, 0.6, 0.6 };
    /* the local frame turned: about x by 50 degrees, then about y by 25 */
    const double ax = 50 * 3.14159265358979 / 180, ay = 25 * 3.14159265358979 / 180;
    const double Rx[3][3] = { { 1, 0, 0 }, { 0, cos(ax), -sin(ax) }, { 0, sin(ax), cos(ax) } };
    const double Ry[3][3] = { { cos(ay), 0, sin(ay) }, { 0, 1, 0 }, { -sin(ay), 0, cos(ay) } };
    bar_scene_t sc = { .t = { -0.4, -0.2, -11.0 }, .radius = radius };
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) { sc.R[i][j] = 0; for (int k = 0; k < 3; ++k) sc.R[i][j] += Ry[i][k] * Rx[k][j]; }
    float pos[N * 3], axis[N * 4], profile[N * 4];
    uint32_t index[N];
    double local[N][3];
    for (int k = 0; k < NA; ++k) { local[k][0] = -3.6 + k * d; local[k][1] = -0.8; local[k][2] = 0; }
    for (int k = 0; k < NB; ++k) { local[NA + k][0] = -2.0 + k * d; local[NA + k][1] = 1.6; local[NA + k][2] = 0; }
    for (int k = 0; k < N; ++k) {
        for (int i = 0; i < 3; ++i) pos[3 * k + i] = (float)(sc.R[i][0] * local[k][0] + sc.R[i][1] * local[k][1] + sc.R[i][2] * local[k][2] + sc.t[i]);
        for (int i = 0; i < 3; ++i) axis[4 * k + i] = (float)(sc.R[i][1] * ((k % 2) ? 2.0 : 0.5));    /* local y, any length */
        axis[4 * k + 3] = 0;
        profile[4 * k] = profile[4 * k + 1] = k < NA ? (float)wa[k] : 0.0f;
        profile[4 * k + 2] = profile[4 * k + 3] = (float)radius;
        index[k] = (uint32_t)k;
    }
    index[0] |= MD_RENDER_SPLINE_FIRST; index[NA - 1] |= MD_RENDER_SPLINE_LAST;
    index[NA] |= MD_RENDER_SPLINE_FIRST; index[N - 1] |= MD_RENDER_SPLINE_LAST;
    uint32_t vis = ~0u & ~(1u << 3);                    /* part 3 hidden */
    bar_part_t parts[N];
    int np = 0;
    for (int k = 0; k < N; ++k) {
        if (k == 3) continue;
        const bool first = k == 0 || k == NA, last = k == NA - 1 || k == N - 1;
        parts[np].x0 = local[k][0] - (first ? 0.0 : 0.5 * d);
        parts[np].x1 = local[k][0] + (last ? 0.0 : 0.5 * d);
        parts[np].yc = local[k][1];
        parts[np].a  = profile[4 * k];
        parts[np].element = k;
        np++;
    }
    sc.part = parts; sc.part_count = np;

    md_render_batch_t b = {0};
    b.key = 9; b.kind = MD_RENDER_SPLINES; b.count = N; b.points.count = N;
    b.points.position = f_buf(&fx, pos, sizeof(pos));
    b.points.visible[0] = f_buf(&fx, &vis, 4);
    b.index = f_buf(&fx, index, sizeof(index));
    b.style.splines.axis    = f_buf(&fx, axis, sizeof(axis));
    b.style.splines.profile = f_buf(&fx, profile, sizeof(profile));
    md_render_frame_t f = frame_of(W, H, cam, &b, 1);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    int judged = 0;
    EXPECT_EQ(0, check_bars(&fx, &c, cam, W, H, &sc, &judged));
    EXPECT_GT(judged, (int)(W * H * 8 / 10));
    printf("    %d parts drawn (one hidden), %d pixels judged\n", np, judged);

    /* the same with the compute raster asked for every part: spline parts
       have none (their code would cost spheres and capsules registers), so
       the same ids and depths to the bit, all from the hardware raster */
    md_render_view_t hw = md_render_view_create(fx.r, "hardware only");
    md_render_frame_t fh = f;
    fh.options = MD_RENDER_NO_COMPUTE_RASTER;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, hw, &fh, NULL));
    f.small_pixels = 1000.0f;
    f.content[MD_RENDER_CONTENT_GEOMETRY]++;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    md_gpu_stream_sync(fx.gfx);
    md_render_stats_t st;
    ASSERT_TRUE(md_render_view_stats(fx.view, &st));
    EXPECT_EQ(0u, st.compute[0] + st.compute[1]);
    EXPECT_EQ((uint32_t)(N - 1), st.hardware[0] + st.hardware[1]);
    int covered = 0;
    EXPECT_EQ(0, compare_views(&fx, fx.view, hw, W, H, &covered));
    md_render_view_destroy(fx.gfx, hw);

    /* a pick at part 5's middle: the part, its point on the face */
    double pl[3] = { local[5][0], local[5][1], radius }, pv[3];
    for (int i = 0; i < 3; ++i) pv[i] = sc.R[i][0] * pl[0] + sc.R[i][1] * pl[1] + sc.R[i][2] * pl[2] + sc.t[i];
    uint32_t qx = (uint32_t)((pv[0] / -pv[2] / (c.tan_half * c.aspect) * 0.5 + 0.5) * W);
    uint32_t qy = (uint32_t)((0.5 - pv[1] / -pv[2] / c.tan_half * 0.5) * H);
    md_render_pick_t pk;
    ASSERT_TRUE(f_pick(&fx, qx, qy, &pk));
    EXPECT_EQ(9u, (unsigned)pk.key);
    EXPECT_EQ(MD_RENDER_SPLINES, pk.kind);
    EXPECT_EQ(5u, pk.element);
    ray_t r = pixel_ray(&c, W, H, qx + 0.5, qy + 0.5);
    double tt = 0.0, gap;
    int who;
    for (int it = 0; it < 4000; ++it) {
        double p[3] = { r.o[0] + r.d[0] * tt, r.o[1] + r.d[1] * tt, r.o[2] + r.d[2] * tt };
        double dd = bar_sdf(&sc, p, &who, &gap);
        if (dd < 1e-7) break;
        tt += dd;
    }
    for (int k = 0; k < 3; ++k) EXPECT_NEAR(r.o[k] + r.d[k] * tt, pk.position[k], 1e-3);
    f_close(&fx);
}


/* Curved splines against the sweep they stand for: the union of the balls
   of the radius around the profile's segments, sampled densely along each
   part from the same B-spline (mirrored ends, unit axes in the sign given,
   made orthogonal to the tangent). The renderer cuts parts into straight
   pieces, so it is judged within a pixel of that surface. */
enum { SW_S = 49 };
typedef struct sweep_part_t {
    int    element;
    double a[SW_S][3], b[SW_S][3];   /* the profile's segment at each sample, view space */
    double r[SW_S];                  /* and its radius */
    double ctr[3], bound;            /* a ball around every segment of the part, radius included */
} sweep_part_t;
typedef struct sweep_t { sweep_part_t* part; int count; } sweep_t;

static double pt_seg_dist(const double p[3], const double a[3], const double b[3]) {
    double ab[3] = { b[0] - a[0], b[1] - a[1], b[2] - a[2] }, ap[3] = { p[0] - a[0], p[1] - a[1], p[2] - a[2] };
    double ll = dot3(ab, ab), s = ll > 0 ? dot3(ap, ab) / ll : 0.0;
    s = s < 0 ? 0 : s > 1 ? 1 : s;
    double d[3] = { ap[0] - ab[0] * s, ap[1] - ab[1] * s, ap[2] - ab[2] * s };
    return sqrt(dot3(d, d));
}

/* Distance between the segments p1 q1 and p2 q2 (Ericson 5.1.9) */
static double seg_seg_dist(const double p1[3], const double q1[3], const double p2[3], const double q2[3]) {
    double d1[3], d2[3], r[3];
    for (int i = 0; i < 3; ++i) { d1[i] = q1[i] - p1[i]; d2[i] = q2[i] - p2[i]; r[i] = p1[i] - p2[i]; }
    double a = dot3(d1, d1), e = dot3(d2, d2), f = dot3(d2, r), s, t;
    if (e <= 1e-30) {
        if (a <= 1e-30) { s = t = 0; }
        else { t = 0; s = -dot3(d1, r) / a; s = s < 0 ? 0 : s > 1 ? 1 : s; }
    } else {
        double c = dot3(d1, r);
        if (a <= 1e-30) { s = 0; t = f / e; t = t < 0 ? 0 : t > 1 ? 1 : t; }
        else {
            double b = dot3(d1, d2), den = a * e - b * b;
            s = den > 1e-30 ? (b * f - c * e) / den : 0.0;
            s = s < 0 ? 0 : s > 1 ? 1 : s;
            t = (b * s + f) / e;
            if (t < 0) { t = 0; s = -c / a; s = s < 0 ? 0 : s > 1 ? 1 : s; }
            else if (t > 1) { t = 1; s = (b - c) / a; s = s < 0 ? 0 : s > 1 ? 1 : s; }
        }
    }
    double dd[3];
    for (int i = 0; i < 3; ++i) dd[i] = p1[i] + d1[i] * s - (p2[i] + d2[i] * t);
    return sqrt(dot3(dd, dd));
}

/* The sweep of every visible part of the splines in pos (view space) */
static int sweep_build(sweep_t* sw, const double (*pos)[3], const double (*axis)[3], const float* profile, const uint32_t* index,
                       const bool* visible, int n) {
    sw->part = calloc(n, sizeof(sweep_part_t));
    sw->count = 0;
    for (int s0 = 0; s0 < n;) {
        int s1 = s0;
        while (!(index[s1] & MD_RENDER_SPLINE_LAST)) s1++;
        int m = s1 - s0 + 1;
        /* the points with two mirrored at each end; axes unit, in the sign given */
        double (*P)[3] = malloc((m + 4) * sizeof(*P)), (*Wx)[3] = malloc((m + 4) * sizeof(*Wx));
        for (int k = 0; k < m; ++k) {
            for (int i = 0; i < 3; ++i) { P[k + 2][i] = pos[s0 + k][i]; Wx[k + 2][i] = axis[s0 + k][i]; }
            double l = sqrt(dot3(Wx[k + 2], Wx[k + 2]));
            for (int i = 0; i < 3; ++i) Wx[k + 2][i] = l > 0 ? Wx[k + 2][i] / l : 0.0;
        }
        for (int i = 0; i < 3; ++i) {
            P[1][i] = 2 * P[2][i] - P[3][i];           P[0][i] = 2 * P[1][i] - P[2][i];
            P[m + 2][i] = 2 * P[m + 1][i] - P[m][i];   P[m + 3][i] = 2 * P[m + 2][i] - P[m + 1][i];
            Wx[0][i] = Wx[1][i] = Wx[2][i];            Wx[m + 2][i] = Wx[m + 3][i] = Wx[m + 1][i];
        }
        for (int k = 0; k < m && m > 1; ++k) {
            int e = s0 + k;
            if (!visible[e]) continue;
            sweep_part_t* pt = &sw->part[sw->count++];
            pt->element = e;
            double u0 = k == 0 ? 0.0 : -0.5, u1 = k == m - 1 ? 0.0 : 0.5;
            for (int j = 0; j < SW_S; ++j) {
                double f = (double)j / (SW_S - 1), u = u0 + (u1 - u0) * f;
                int base = u < 0 ? k : k + 1;          /* extended index of the first of the four */
                double x = u < 0 ? u + 1.0 : u, x2 = x * x, x3 = x2 * x, ix = 1 - x;
                double bw[4] = { ix * ix * ix / 6, (3 * x3 - 6 * x2 + 4) / 6, (-3 * x3 + 3 * x2 + 3 * x + 1) / 6, x3 / 6 };
                double dw[4] = { -0.5 * ix * ix, 1.5 * x2 - 2 * x, -1.5 * x2 + x + 0.5, 0.5 * x2 };
                double c[3] = {0}, t[3] = {0}, w[3] = {0};
                for (int q = 0; q < 4; ++q) for (int i = 0; i < 3; ++i) {
                    c[i] += P[base + q][i] * bw[q]; t[i] += P[base + q][i] * dw[q]; w[i] += Wx[base + q][i] * bw[q];
                }
                double tt = dot3(t, t), wt = dot3(w, t);
                if (tt > 0) for (int i = 0; i < 3; ++i) w[i] -= t[i] * wt / tt;
                double wl = sqrt(dot3(w, w));
                double a = (double)profile[4 * e] + ((double)profile[4 * e + 1] - (double)profile[4 * e]) * f;
                pt->r[j] = (double)profile[4 * e + 2] + ((double)profile[4 * e + 3] - (double)profile[4 * e + 2]) * f;
                for (int i = 0; i < 3; ++i) {
                    double wi = wl > 1e-12 ? w[i] / wl : 0.0;
                    pt->a[j][i] = c[i] - wi * a;
                    pt->b[j][i] = c[i] + wi * a;
                }
            }
            for (int i = 0; i < 3; ++i) pt->ctr[i] = 0.5 * (pt->a[SW_S / 2][i] + pt->b[SW_S / 2][i]);
            pt->bound = 0;
            for (int j = 0; j < SW_S; ++j) {
                double da[3] = { pt->a[j][0] - pt->ctr[0], pt->a[j][1] - pt->ctr[1], pt->a[j][2] - pt->ctr[2] };
                double db[3] = { pt->b[j][0] - pt->ctr[0], pt->b[j][1] - pt->ctr[1], pt->b[j][2] - pt->ctr[2] };
                pt->bound = fmax(pt->bound, fmax(sqrt(dot3(da, da)), sqrt(dot3(db, db))) + pt->r[j]);
            }
        }
        free(P); free(Wx);
        s0 = s1 + 1;
    }
    return sw->count;
}

/* Distance from p to the sweep; the nearest part's element and how much farther the next is */
static double sweep_sdf(const sweep_t* sw, const double p[3], int* who, double* gap) {
    double best = 1e30, second = 1e30;
    *who = -1;
    for (int k = 0; k < sw->count; ++k) {
        const sweep_part_t* pt = &sw->part[k];
        double dc[3] = { p[0] - pt->ctr[0], p[1] - pt->ctr[1], p[2] - pt->ctr[2] };
        if (sqrt(dot3(dc, dc)) - pt->bound > second) continue;
        double d = 1e30;
        for (int j = 0; j < SW_S; ++j) d = fmin(d, pt_seg_dist(p, pt->a[j], pt->b[j]) - pt->r[j]);
        if (d < best) { second = best; best = d; *who = pt->element; }
        else if (d < second) second = d;
    }
    *gap = second - best;
    return best;
}

/* The sweep's normal at p (on its surface) by part k's nearest point:
   between the samples nearest it the segments are interpolated, and the
   nearest of them found by ternary search, good to ~0.1 degree here */
static void sweep_normal(const sweep_t* sw, int k, const double p[3], double n[3]) {
    const sweep_part_t* pt = &sw->part[k];
    int best = 0;
    double bd = 1e30;
    for (int j = 0; j < SW_S; ++j) {
        double d = pt_seg_dist(p, pt->a[j], pt->b[j]) - pt->r[j];
        if (d < bd) { bd = d; best = j; }
    }
    double lo = best > 0 ? best - 1 : 0, hi = best < SW_S - 1 ? best + 1 : SW_S - 1;
    double sa[3], sb[3], sr;
    #define SEG_AT(x) do { int j0 = (int)floor(x); if (j0 >= SW_S - 1) j0 = SW_S - 2; double f = (x) - j0; \
        for (int i = 0; i < 3; ++i) { sa[i] = pt->a[j0][i] + (pt->a[j0 + 1][i] - pt->a[j0][i]) * f; sb[i] = pt->b[j0][i] + (pt->b[j0 + 1][i] - pt->b[j0][i]) * f; } \
        sr = pt->r[j0] + (pt->r[j0 + 1] - pt->r[j0]) * f; } while (0)
    for (int it = 0; it < 60; ++it) {
        double m1 = lo + (hi - lo) / 3, m2 = hi - (hi - lo) / 3;
        SEG_AT(m1); double d1 = pt_seg_dist(p, sa, sb) - sr;
        SEG_AT(m2); double d2 = pt_seg_dist(p, sa, sb) - sr;
        if (d1 < d2) hi = m2; else lo = m1;
    }
    SEG_AT(0.5 * (lo + hi));
    #undef SEG_AT
    double ab[3] = { sb[0] - sa[0], sb[1] - sa[1], sb[2] - sa[2] }, ap[3] = { p[0] - sa[0], p[1] - sa[1], p[2] - sa[2] };
    double ll = dot3(ab, ab), t = ll > 0 ? dot3(ap, ab) / ll : 0.0;
    t = t < 0 ? 0 : t > 1 ? 1 : t;
    for (int i = 0; i < 3; ++i) n[i] = ap[i] - ab[i] * t;
    double l = sqrt(dot3(n, n));
    for (int i = 0; i < 3; ++i) n[i] /= l;
}

static int cmp_double(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return x < y ? -1 : x > y ? 1 : 0;
}

static int sweep_part_of(const sweep_t* sw, int element) {
    for (int k = 0; k < sw->count; ++k) if (sw->part[k].element == element) return k;
    return -1;
}

static float half_bits_to_float(uint32_t h) {
    uint32_t s = (h >> 15) & 1u, e = (h >> 10) & 31u, m = h & 1023u;
    float v = e == 0 ? ldexpf((float)m, -24) : e == 31 ? INFINITY : ldexpf((float)(m | 1024u), (int)e - 25);
    return s ? -v : v;
}

static void oct_decode(uint32_t bits, double n[3]) {
    double x = half_bits_to_float(bits & 0xFFFFu), y = half_bits_to_float(bits >> 16), z = 1.0 - fabs(x) - fabs(y);
    if (z < 0) { double px = (1 - fabs(y)) * (x >= 0 ? 1 : -1), py = (1 - fabs(x)) * (y >= 0 ? 1 : -1); x = px; y = py; }
    double l = sqrt(x * x + y * y + z * z);
    n[0] = x / l; n[1] = y / l; n[2] = z / l;
}

/* How far the ray (t in [0, tmax]) passes outside the sweep: < 0, through it */
static double sweep_ray_clearance(const sweep_t* sw, ray_t r, double tmax) {
    double q[3] = { r.o[0] + r.d[0] * tmax, r.o[1] + r.d[1] * tmax, r.o[2] + r.d[2] * tmax };
    double best = 1e30;
    for (int k = 0; k < sw->count; ++k) {
        const sweep_part_t* pt = &sw->part[k];
        if (pt_seg_dist(pt->ctr, r.o, q) - pt->bound > best) continue;
        for (int j = 0; j < SW_S; ++j) best = fmin(best, seg_seg_dist(r.o, q, pt->a[j], pt->b[j]) - pt->r[j]);
    }
    return best;
}

UTEST(render, spline_curves_match_the_sweep) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 192, H = 144;
    cam_t c = { .tan_half = 0.45f, .aspect = (float)W / H, .n = 0.1f };
    md_render_camera_t cam = make_camera(&c);
    /* A: a helix ribbon, axes along the helix axis in alternating length, out of a coil over its first two parts (a tube tapering into
       the ribbon); B: a twisting strand with an arrow (a step up, then down
       to 0, thinner); C: a wavy tube, thicker along it, one part hidden */
    enum { NA = 11, NB = 8, NC = 9, N = NA + NB + NC };
    double local[N][3], laxis[N][3];
    float profile[N * 4];
    #define PROFILE(e, a0, a1, r0, r1) (profile[4 * (e)] = (float)(a0), profile[4 * (e) + 1] = (float)(a1), \
                                        profile[4 * (e) + 2] = (float)(r0), profile[4 * (e) + 3] = (float)(r1))
    for (int k = 0; k < NA; ++k) {
        double ph = k * 100.0 * 3.14159265358979 / 180;
        local[k][0] = -5.0 + k * 1.0; local[k][1] = 2.6 + 2.0 * cos(ph); local[k][2] = 2.0 * sin(ph);
        double l = (k % 2) ? 0.5 : 1.5;
        laxis[k][0] = l; laxis[k][1] = 0; laxis[k][2] = 0;
        PROFILE(k, 0.7, 0.7, 0.18, 0.18);
    }
    PROFILE(0, 0.0, 0.35, 0.3, 0.24);
    PROFILE(1, 0.35, 0.7, 0.24, 0.18);
    for (int k = 0; k < NB; ++k) {
        int e = NA + k;
        double th = k * 15.0 * 3.14159265358979 / 180, l = 1.0 + k % 3;
        local[e][0] = -5.0 + k * 1.4; local[e][1] = -1.4 + ((k % 2) ? 0.35 : -0.35); local[e][2] = 0.3 * sin((double)k);
        laxis[e][0] = 0; laxis[e][1] = sin(th) * l; laxis[e][2] = cos(th) * l;
        PROFILE(e, 0.8, 0.8, 0.18, 0.18);
    }
    PROFILE(NA + NB - 2, 1.3, 0.65, 0.24, 0.2);                 /* the arrow */
    PROFILE(NA + NB - 1, 0.65, 0.0, 0.2, 0.16);
    for (int k = 0; k < NC; ++k) {
        int e = NA + NB + k;
        local[e][0] = -5.0 + k * 1.25; local[e][1] = -4.3 + 0.8 * sin(1.3 * k); local[e][2] = 1.0 * cos(0.9 * k);
        laxis[e][0] = laxis[e][1] = laxis[e][2] = 0;
        double x0 = k == 0 ? 0.0 : k - 0.5, x1 = k == NC - 1 ? k : k + 0.5;
        PROFILE(e, 0.0, 0.0, 0.1 + 0.2 * x0 / (NC - 1), 0.1 + 0.2 * x1 / (NC - 1));
    }
    #undef PROFILE
    const double ax = 35 * 3.14159265358979 / 180, ay = -20 * 3.14159265358979 / 180;
    const double Rx[3][3] = { { 1, 0, 0 }, { 0, cos(ax), -sin(ax) }, { 0, sin(ax), cos(ax) } };
    const double Ry[3][3] = { { cos(ay), 0, sin(ay) }, { 0, 1, 0 }, { -sin(ay), 0, cos(ay) } };
    const double tr[3] = { 0.3, 0.2, -15.0 };
    double R[3][3];
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) { R[i][j] = 0; for (int k = 0; k < 3; ++k) R[i][j] += Ry[i][k] * Rx[k][j]; }
    float pos[N * 3], axis[N * 4];
    double vpos[N][3], vaxis[N][3];
    uint32_t index[N];
    bool visible[N];
    for (int k = 0; k < N; ++k) {
        for (int i = 0; i < 3; ++i) {
            pos[3 * k + i] = (float)(R[i][0] * local[k][0] + R[i][1] * local[k][1] + R[i][2] * local[k][2] + tr[i]);
            axis[4 * k + i] = (float)(R[i][0] * laxis[k][0] + R[i][1] * laxis[k][1] + R[i][2] * laxis[k][2]);
            vpos[k][i] = pos[3 * k + i];
            vaxis[k][i] = axis[4 * k + i];
        }
        axis[4 * k + 3] = 0;
        index[k] = (uint32_t)k;
        visible[k] = true;
    }
    index[0] |= MD_RENDER_SPLINE_FIRST;       index[NA - 1] |= MD_RENDER_SPLINE_LAST;
    index[NA] |= MD_RENDER_SPLINE_FIRST;      index[NA + NB - 1] |= MD_RENDER_SPLINE_LAST;
    index[NA + NB] |= MD_RENDER_SPLINE_FIRST; index[N - 1] |= MD_RENDER_SPLINE_LAST;
    const int hidden = NA + NB + 4;
    visible[hidden] = false;
    uint32_t vis = ~(1u << hidden);

    md_render_batch_t b = {0};
    b.key = 5; b.kind = MD_RENDER_SPLINES; b.count = N; b.points.count = N;
    b.points.position = f_buf(&fx, pos, sizeof(pos));
    b.points.visible[0] = f_buf(&fx, &vis, 4);
    b.index = f_buf(&fx, index, sizeof(index));
    b.style.splines.axis    = f_buf(&fx, axis, sizeof(axis));
    b.style.splines.profile = f_buf(&fx, profile, sizeof(profile));
    md_render_frame_t f = frame_of(W, H, cam, &b, 1);
    f.options = MD_RENDER_NO_COMPUTE_RASTER;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));

    sweep_t sw;
    sweep_build(&sw, (const double (*)[3])vpos, (const double (*)[3])vaxis, profile, index, visible, N);
    md_render_targets_t tg = md_render_view_targets(fx.view);
    uint32_t* id = malloc(W * H * 8);
    float*    dp = malloc(W * H * 4);
    ASSERT_TRUE(f_read(&fx, tg.id, id, W * H * 8) && f_read(&fx, tg.depth, dp, W * H * 4));
    int bad = 0, covered = 0, judged_parts = 0;
    double worst_surface = 0;    /* in pixels */
    for (uint32_t y = 0; y < H; ++y) {
        for (uint32_t x = 0; x < W; ++x) {
            ray_t r = pixel_ray(&c, W, H, x + 0.5, y + 0.5);
            uint32_t ix = id[2 * (y * W + x)], iy = id[2 * (y * W + x) + 1];
            double clr = sweep_ray_clearance(&sw, r, 100.0);
            /* a pixel's size in world at depth z, from the ray through the scene's middle */
            double px_at15 = 15.0 * 2.0 * c.tan_half / H;
            if (iy == 0) {
                if (clr < -0.75 * px_at15) { bad++; if (bad < 8) printf("    pixel (%u, %u): empty, the sweep %.3f px inside\n", x, y, -clr / px_at15); }
                continue;
            }
            covered++;
            if (clr > 0.75 * px_at15) { bad++; if (bad < 8) printf("    pixel (%u, %u): covered, the sweep %.3f px away\n", x, y, clr / px_at15); continue; }
            double zv = -c.n / dp[y * W + x], tp = zv / r.d[2];
            double P[3] = { r.o[0] + r.d[0] * tp, r.o[1] + r.d[1] * tp, r.o[2] + r.d[2] * tp };
            double pxw = -P[2] * 2.0 * c.tan_half / H;
            int who; double gap;
            double sd = sweep_sdf(&sw, P, &who, &gap);
            worst_surface = fmax(worst_surface, fabs(sd) / pxw);
            if (fabs(sd) > 0.5 * pxw) { bad++; if (bad < 8) printf("    pixel (%u, %u): %.3f px off the surface\n", x, y, sd / pxw); continue; }
            if (gap > pxw) {
                judged_parts++;
                if (ix != (uint32_t)who) { bad++; if (bad < 8) printf("    pixel (%u, %u): part %u, expected %d\n", x, y, ix, who); continue; }
            }
            /* nothing of the sweep in front: the ray before P stays outside but for a pixel */
            double tt = 0.0;
            for (int it = 0; it < 400 && tt < tp; ++it) {
                double q[3] = { r.o[0] + r.d[0] * tt, r.o[1] + r.d[1] * tt, r.o[2] + r.d[2] * tt };
                double d = sweep_sdf(&sw, q, &who, &gap);
                double qw = -q[2] * 2.0 * c.tan_half / H;
                if (d < -0.5 * qw) { bad++; if (bad < 8) printf("    pixel (%u, %u): the sweep %.3f px in front\n", x, y, (tp - tt) / qw); break; }
                tt += fmax(fabs(d), 0.25 * qw);
            }
        }
    }
    printf("    %d pixels covered, %d with a clear part; the farthest %.3f px from the sweep\n", covered, judged_parts, worst_surface);

    /* Normals against the sweep's, where the part is clear: smooth across
       the faces' edges into the rims, and along the pieces. The faces were
       once built on the centre's normal and met twisting rims at a crease
       (13 degrees on this helix, a groove under SSAO): its pixels stand
       out at the top of this distribution. */
    {
        uint32_t* nb = malloc(W * H * 4);
        ASSERT_TRUE(f_read(&fx, tg.normal, nb, W * H * 4));
        double* errs = malloc(W * H * sizeof(double));
        int ne = 0;
        for (uint32_t y = 0; y < H; ++y) {
            for (uint32_t x = 0; x < W; ++x) {
                uint32_t i = y * W + x;
                if (id[2 * i + 1] == 0) continue;
                ray_t r = pixel_ray(&c, W, H, x + 0.5, y + 0.5);
                double zv = -c.n / dp[i], tp = zv / r.d[2];
                double P[3] = { r.o[0] + r.d[0] * tp, r.o[1] + r.d[1] * tp, r.o[2] + r.d[2] * tp };
                double pxw = -P[2] * 2.0 * c.tan_half / H;
                int who; double gap;
                sweep_sdf(&sw, P, &who, &gap);
                if (gap < 2.0 * pxw || who != (int)id[2 * i]) continue;
                double nr[3], ng[3];
                sweep_normal(&sw, sweep_part_of(&sw, who), P, nr);
                if (fabs(dot3(nr, r.d)) < 0.3) continue;   /* grazing: the depth's rounding moves P around a thin rim */
                oct_decode(nb[i], ng);
                double cs = dot3(nr, ng);
                errs[ne++] = acos(cs > 1 ? 1 : cs < -1 ? -1 : cs) * 180.0 / 3.14159265358979;
            }
        }
        qsort(errs, (size_t)ne, sizeof(double), cmp_double);
        double p50 = ne ? errs[ne / 2] : 0, p99 = ne ? errs[(ne * 99) / 100] : 0, pmax = ne ? errs[ne - 1] : 0;
        printf("    normals against the sweep's, %d pixels: median %.2f, 99th %.2f, largest %.2f degrees\n", ne, p50, p99, pmax);
        EXPECT_GT(ne, judged_parts / 2);
        /* The rest is the pieces' sag over a thin rim's radius, near where
           it meets the face (the face's normal is the curve's, the rim's
           point is the piece's); it shrinks as the view comes closer. The
           faces on the centre's normal: median 4.8, 99th 30 degrees. */
        EXPECT_LT(p50, 1.0);
        EXPECT_LT(p99, 6.0);
        free(errs); free(nb);
    }
    EXPECT_EQ(0, bad);
    EXPECT_GT(covered, (int)(W * H / 8));
    EXPECT_GT(judged_parts, covered / 2);

    /* the compute raster asked for: spline parts stay with the hardware
       raster, the same ids and depths to the bit */
    md_render_view_t cr = md_render_view_create(fx.r, "compute raster");
    md_render_frame_t fc = f;
    fc.options = 0;
    fc.small_pixels = 1000.0f;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, cr, &fc, NULL));
    md_gpu_stream_sync(fx.gfx);
    md_render_stats_t st;
    ASSERT_TRUE(md_render_view_stats(cr, &st));
    EXPECT_EQ(0u, st.compute[0] + st.compute[1]);
    EXPECT_GT(st.hardware[0] + st.hardware[1], (uint32_t)(N / 2));
    int cov2 = 0;
    EXPECT_EQ(0, compare_views(&fx, cr, fx.view, W, H, &cov2));
    md_render_view_destroy(fx.gfx, cr);

    /* picks at each part's control point: what the raster drew there */
    int picks = 0;
    for (int e = 0; e < N; ++e) {
        double* q = vpos[e];
        int32_t qx = (int32_t)((q[0] / -q[2] / (c.tan_half * c.aspect) * 0.5 + 0.5) * W);
        int32_t qy = (int32_t)((0.5 - q[1] / -q[2] / c.tan_half * 0.5) * H);
        if (qx < 0 || qy < 0 || qx >= (int32_t)W || qy >= (int32_t)H) continue;
        uint32_t i = (uint32_t)qy * W + (uint32_t)qx;
        md_render_pick_t pk;
        ASSERT_TRUE(f_pick(&fx, (uint32_t)qx, (uint32_t)qy, &pk));
        if (id[2 * i + 1] == 0) { EXPECT_EQ(0u, (unsigned)pk.key); continue; }
        picks++;
        EXPECT_EQ(5u, (unsigned)pk.key);
        EXPECT_EQ(id[2 * i], pk.element);
        EXPECT_EQ(pk.element & 0x3FFFFFFFu, pk.attr);
        ray_t r = pixel_ray(&c, W, H, qx + 0.5, qy + 0.5);
        double zv = -c.n / dp[i], tp = zv / r.d[2];
        for (int k = 0; k < 3; ++k) EXPECT_NEAR(r.o[k] + r.d[k] * tp, pk.position[k], 1e-4 * 15.0);
    }
    EXPECT_GT(picks, N / 2);

    /* every other part of the helix selected: its mask where the raster
       shows them, its outline around more (where the turn in front covers
       the turn behind) */
    uint32_t sel = 0x555u & ((1u << NA) - 1u);
    b.points.state[0] = f_buf(&fx, &sel, 4);
    f.content[MD_RENDER_CONTENT_STATE]++;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    uint32_t* m = malloc(W * H * 4);
    ASSERT_TRUE(f_read(&fx, tg.mask, m, W * H * 4));
    int mbad = 0, mvis = 0, many = 0;
    for (uint32_t i = 0; i < W * H; ++i) {
        bool helix = id[2 * i + 1] != 0 && id[2 * i] < NA && (id[2 * i] % 2) == 0;
        bool vis = (m[i] & 0xFFu) != 0, any = (m[i] & 0xFF00u) != 0;
        mvis += vis; many += any;
        if (vis != helix || (vis && !any) || (m[i] & 0xFFFF0000u)) {
            if (mbad < 6) printf("    pixel (%u, %u): mask %08x, part %u\n", i % W, i / W, m[i], id[2 * i]);
            mbad++;
        }
    }
    EXPECT_EQ(0, mbad);
    EXPECT_GT(mvis, 500);
    EXPECT_GT(many, mvis + 50);
    printf("    the selected parts: %d pixels in view, %d in all\n", mvis, many);
    free(m);

    /* shrunk across a fade band through the scene: thinner, the rasters
       still agree to the bit; past the band, nothing */
    b.points.state[0] = 0;
    b.fade.mode = MD_RENDER_FADE_SHRINK;
    b.fade.out[0] = 14.0f; b.fade.out[1] = 8.0f;     /* px per unit: depth 11.4 to 20 */
    f.content[MD_RENDER_CONTENT_GEOMETRY]++;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    cr = md_render_view_create(fx.r, "compute raster");
    ASSERT_TRUE(md_render_view_draw(fx.gfx, cr, &fc, NULL));
    md_gpu_stream_sync(fx.gfx);
    int cov3 = 0;
    EXPECT_EQ(0, compare_views(&fx, cr, fx.view, W, H, &cov3));
    EXPECT_GT(cov3, covered / 4);
    EXPECT_LT(cov3, covered * 3 / 4);
    b.fade.out[0] = 40.0f; b.fade.out[1] = 30.0f;    /* depth 4 to 5.3: the scene is past it */
    f.content[MD_RENDER_CONTENT_GEOMETRY]++;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    md_gpu_stream_sync(fx.gfx);
    ASSERT_TRUE(md_render_view_stats(fx.view, &st));
    EXPECT_EQ((uint32_t)(N - 1), st.faded);
    ASSERT_TRUE(f_read(&fx, tg.id, id, W * H * 8));
    int left = 0;
    for (uint32_t i = 0; i < W * H; ++i) left += id[2 * i + 1] != 0;
    EXPECT_EQ(0, left);
    printf("    shrunk: %d pixels; past the band: %d\n", cov3, left);
    md_render_view_destroy(fx.gfx, cr);

    /* with spheres in the same frame: the main kernels and the spline twins
       each take theirs, nothing twice, nothing lost. Every pixel as the
       nearer of the two drawn alone: depth, id, colour and normal. */
    b.fade.mode = 0;
    enum { NS = 24 };
    float spos[NS * 3];
    for (int k = 0; k < NS; ++k) {
        int e = (k * 7) % N;                                /* by control points, in front and behind */
        for (int i = 0; i < 3; ++i) spos[3 * k + i] = pos[3 * e + i];
        spos[3 * k + 0] += 0.3f * (float)((k % 3) - 1);
        spos[3 * k + 2] += (k % 2) ? 0.6f : -0.6f;
    }
    md_render_batch_t bs = {0};
    bs.key = 6; bs.kind = MD_RENDER_SPHERES; bs.count = NS; bs.points.count = NS;
    bs.points.position = f_buf(&fx, spos, sizeof(spos));
    bs.style.spheres.radius_add = 0.45f;
    bs.color.mode = MD_RENDER_COLOR_UNIFORM;
    bs.color.uniform = 0xFF30C030u;
    md_render_batch_t both[2] = { bs, b };
    md_render_view_t vs = md_render_view_create(fx.r, "spheres"), vp = md_render_view_create(fx.r, "splines"), vb = md_render_view_create(fx.r, "both");
    md_render_frame_t fs = frame_of(W, H, cam, &bs, 1), fp = frame_of(W, H, cam, &b, 1), fbo = frame_of(W, H, cam, both, 2);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, vs, &fs, NULL));
    ASSERT_TRUE(md_render_view_draw(fx.gfx, vp, &fp, NULL));
    ASSERT_TRUE(md_render_view_draw(fx.gfx, vb, &fbo, NULL));
    md_gpu_stream_sync(fx.gfx);
    md_render_view_t views[3] = { vs, vp, vb };
    uint32_t* ids[3]; float* dps[3]; uint32_t* cols[3]; uint32_t* nrms[3];
    for (int k = 0; k < 3; ++k) {
        md_render_targets_t t = md_render_view_targets(views[k]);
        ids[k] = malloc(W * H * 8); dps[k] = malloc(W * H * 4); cols[k] = malloc(W * H * 4); nrms[k] = malloc(W * H * 4);
        ASSERT_TRUE(f_read(&fx, t.id, ids[k], W * H * 8) && f_read(&fx, t.depth, dps[k], W * H * 4) &&
                    f_read(&fx, t.color, cols[k], W * H * 4) && f_read(&fx, t.normal, nrms[k], W * H * 4));
    }
    int mixed_bad = 0, from[2] = { 0, 0 };
    for (uint32_t i = 0; i < W * H; ++i) {
        bool hs = ids[0][2 * i + 1] != 0, hp = ids[1][2 * i + 1] != 0;
        if (!hs && !hp) { if (ids[2][2 * i + 1] != 0) mixed_bad++; continue; }
        int w = !hp || (hs && dps[0][i] >= dps[1][i]) ? 0 : 1;      /* reverse Z: the larger is nearer */
        from[w]++;
        uint32_t ey = (ids[w][2 * i + 1] & 0xFFFFFu) | ((uint32_t)(w + 1) << 20);   /* slot 0 alone, slot w in both */
        if (ids[2][2 * i] != ids[w][2 * i] || ids[2][2 * i + 1] != ey || dps[2][i] != dps[w][i] ||
            cols[2][i] != cols[w][i] || nrms[2][i] != nrms[w][i]) {
            if (mixed_bad < 6) printf("    pixel (%u, %u): both %u/%x, alone (%s) %u/%x\n", i % W, i / W, ids[2][2 * i], ids[2][2 * i + 1],
                                      w ? "splines" : "spheres", ids[w][2 * i], ids[w][2 * i + 1]);
            mixed_bad++;
        }
    }
    EXPECT_EQ(0, mixed_bad);
    EXPECT_GT(from[0], 200);
    EXPECT_GT(from[1], 1000);
    printf("    with spheres: %d pixels theirs, %d the splines'\n", from[0], from[1]);
    for (int k = 0; k < 3; ++k) { free(ids[k]); free(dps[k]); free(cols[k]); free(nrms[k]); md_render_view_destroy(fx.gfx, views[k]); }
    free(id); free(dp); free(sw.part);
    f_close(&fx);
}

/* ---- Polygons ------------------------------------------------------------------ */

/* A polygon as md_render.h defines it, in double: the corners' centroid,
   Newell's normal, the corners projected onto that plane; rounded by r. On a
   line: its edges alone. */
typedef struct poly_ref_t {
    double q[8][3], c[3], n[3], r, bound;
    int    corners;
    bool   plate;
    uint32_t key, linear;                /* its batch, and its id's element (instance * count + element) */
    uint32_t point[8];                   /* its corners' points: their attributes */
} poly_ref_t;

static void cross3(const double* a, const double* b, double* o) {
    o[0] = a[1] * b[2] - a[2] * b[1]; o[1] = a[2] * b[0] - a[0] * b[2]; o[2] = a[0] * b[1] - a[1] * b[0];
}

static void poly_ref_build(poly_ref_t* p, const double (*q)[3], int n, double r) {
    p->corners = n;
    p->r = r;
    for (int i = 0; i < 3; ++i) { p->c[i] = 0; for (int k = 0; k < n; ++k) p->c[i] += q[k][i] / n; }
    double nv[3] = { 0, 0, 0 }, sz = 0;
    for (int k = 0; k < n; ++k) {
        const double* b = q[(k + 1) % n];
        double u[3] = { q[k][0] - p->c[0], q[k][1] - p->c[1], q[k][2] - p->c[2] }, v[3] = { b[0] - p->c[0], b[1] - p->c[1], b[2] - p->c[2] }, x[3];
        cross3(u, v, x);
        for (int i = 0; i < 3; ++i) nv[i] += x[i];
        sz = fmax(sz, sqrt(dot3(u, u)));
    }
    double nl = sqrt(dot3(nv, nv));
    p->plate = nl > 1e-6 * sz * sz;
    for (int i = 0; i < 3; ++i) p->n[i] = p->plate ? nv[i] / nl : 0.0;
    p->bound = 0;
    for (int k = 0; k < n; ++k) {
        double d[3] = { q[k][0] - p->c[0], q[k][1] - p->c[1], q[k][2] - p->c[2] }, h = dot3(d, p->n);
        for (int i = 0; i < 3; ++i) p->q[k][i] = q[k][i] - p->n[i] * h;
        double e[3] = { p->q[k][0] - p->c[0], p->q[k][1] - p->c[1], p->q[k][2] - p->c[2] };
        p->bound = fmax(p->bound, sqrt(dot3(e, e)));
    }
    p->bound += r;
}

/* Whether x, on the plane, is inside the polygon */
static bool poly_ref_inside(const poly_ref_t* p, const double x[3]) {
    if (!p->plate) return false;
    for (int k = 0; k < p->corners; ++k) {
        const double* a = p->q[k]; const double* b = p->q[(k + 1) % p->corners];
        double e[3] = { b[0] - a[0], b[1] - a[1], b[2] - a[2] }, d[3] = { x[0] - a[0], x[1] - a[1], x[2] - a[2] }, s[3];
        cross3(e, d, s);
        if (dot3(s, p->n) < 0) return false;
    }
    return true;
}

/* Distance from x to the flat polygon (not rounded), and its nearest point */
static double poly_ref_dist(const poly_ref_t* p, const double x[3], double cp[3]) {
    double d[3] = { x[0] - p->c[0], x[1] - p->c[1], x[2] - p->c[2] }, h = dot3(d, p->n);
    double xp[3] = { x[0] - p->n[0] * h, x[1] - p->n[1] * h, x[2] - p->n[2] * h };
    if (poly_ref_inside(p, xp)) { memcpy(cp, xp, sizeof(xp)); return fabs(h); }
    double best = 1e30;
    for (int k = 0; k < p->corners; ++k) {
        const double* a = p->q[k]; const double* b = p->q[(k + 1) % p->corners];
        double ab[3] = { b[0] - a[0], b[1] - a[1], b[2] - a[2] }, ax[3] = { x[0] - a[0], x[1] - a[1], x[2] - a[2] };
        double ll = dot3(ab, ab), s = ll > 0 ? dot3(ax, ab) / ll : 0;
        s = s < 0 ? 0 : s > 1 ? 1 : s;
        double q[3] = { a[0] + ab[0] * s, a[1] + ab[1] * s, a[2] + ab[2] * s }, e[3] = { x[0] - q[0], x[1] - q[1], x[2] - q[2] };
        double dd = sqrt(dot3(e, e));
        if (dd < best) { best = dd; memcpy(cp, q, sizeof(q)); }
    }
    return best;
}

/* The signed distance to the nearest polygon, which, and how much farther the next one is */
static double polys_sdf(const poly_ref_t* P, int count, const double x[3], int* who, double* gap) {
    double best = 1e30, second = 1e30, cp[3];
    *who = -1;
    for (int k = 0; k < count; ++k) {
        double d = poly_ref_dist(&P[k], x, cp) - P[k].r;
        if (d < best) { second = best; best = d; *who = k; }
        else if (d < second) second = d;
    }
    *gap = second - best;
    return best;
}

/* How far the ray (t in [0, tmax]) passes outside every polygon: < 0, through one */
static double polys_ray_clearance(const poly_ref_t* P, int count, ray_t r, double tmax) {
    double q[3] = { r.o[0] + r.d[0] * tmax, r.o[1] + r.d[1] * tmax, r.o[2] + r.d[2] * tmax }, cp[3];
    double best = 1e30;
    for (int k = 0; k < count; ++k) {
        const poly_ref_t* p = &P[k];
        if (pt_seg_dist(p->c, r.o, q) - p->bound > best) continue;
        /* the segment and the convex polygon: through it, or nearest at an edge or an end */
        double d = fmin(poly_ref_dist(p, r.o, cp), poly_ref_dist(p, q, cp));
        double o[3] = { r.o[0] - p->c[0], r.o[1] - p->c[1], r.o[2] - p->c[2] }, e[3] = { q[0] - p->c[0], q[1] - p->c[1], q[2] - p->c[2] };
        double h0 = dot3(o, p->n), h1 = dot3(e, p->n);
        if (p->plate && (h0 > 0) != (h1 > 0)) {
            double s = h0 / (h0 - h1), x[3] = { r.o[0] + (q[0] - r.o[0]) * s, r.o[1] + (q[1] - r.o[1]) * s, r.o[2] + (q[2] - r.o[2]) * s };
            if (poly_ref_inside(p, x)) d = 0;
        }
        for (int j = 0; j < p->corners && d > 0; ++j) d = fmin(d, seg_seg_dist(r.o, q, p->q[j], p->q[(j + 1) % p->corners]));
        best = fmin(best, d - p->r);
    }
    return best;
}

/* A regular polygon of n corners and circumradius rho about ctr, its plane
   tilted by `tilt` degrees from facing the camera about the axis at `spin`
   degrees in the image, wound either way; corners off the plane by +-bump
   in turn */
static void poly_corners(double (*out)[3], int n, const double ctr[3], double rho, double tilt, double spin, bool reverse, double bump) {
    double tl = tilt * 3.14159265358979 / 180, sp = spin * 3.14159265358979 / 180;
    double ax[3] = { cos(sp), sin(sp), 0 };                  /* the tilt's axis, in the image plane */
    double nz[3] = { 0, 0, 1 }, w[3];
    cross3(nz, ax, w);                                       /* in the image plane, across the axis */
    double u[3], v[3], nn[3];
    for (int i = 0; i < 3; ++i) { u[i] = ax[i]; v[i] = w[i] * cos(tl) + nz[i] * sin(tl); }
    cross3(u, v, nn);
    for (int k = 0; k < n; ++k) {
        int j = reverse ? n - 1 - k : k;
        double th = 2 * 3.14159265358979 * j / n + 0.3;
        double z = (k % 2 ? 1.0 : -1.0) * bump;
        for (int i = 0; i < 3; ++i) out[k][i] = ctr[i] + rho * (cos(th) * u[i] + sin(th) * v[i]) + z * nn[i];
    }
}

UTEST(render, polygons_match_their_solids) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 192, H = 144;
    cam_t c = { .tan_half = 0.45f, .aspect = (float)W / H, .n = 0.1f };
    md_render_camera_t cam = make_camera(&c);

    /* Hexagons through an index, one bent out of its plane, one edge-on,
       one with a hidden corner; pentagons in order, twice (an instance
       behind, turned and scaled 1.1); two triangles, one on a line; an
       octagon */
    enum { NH = 6, NP = 4, NT = 2, NO = 1 };
    double hc[NH][3] = { { -6.0, 3.4, -15.0 }, { -2.2, 3.6, -14.0 }, { 1.8, 3.4, -15.5 }, { -6.0, 0.2, -14.5 }, { -2.0, 0.0, -15.0 }, { 2.0, 0.2, -14.0 } };
    double htilt[NH] = { 0, 35, 60, 80, 90, 50 }, hspin[NH] = { 0, 0, 90, 30, 0, 120 };
    double pc[NP][3] = { { -6.5, -3.4, -15.0 }, { -3.0, -3.3, -14.5 }, { 0.5, -3.5, -15.5 }, { 4.0, -3.2, -15.0 } };
    double ptilt[NP] = { 20, 70, 0, 45 }, pspin[NP] = { 60, 10, 0, 200 };
    float  hpos[NH * 6 * 3], ppos[NP * 5 * 3], tpos[NT * 3 * 3], opos[NO * 8 * 3];
    uint32_t hindex[NH * 6];
    double q[8][3];
    for (int k = 0; k < NH; ++k) {
        poly_corners(q, 6, hc[k], 1.3, htilt[k], hspin[k], k == 2, k == 0 ? 0.12 : 0.0);
        for (int j = 0; j < 6; ++j) {
            uint32_t pnt = (uint32_t)(NH * 6 - 1 - (k * 6 + j));    /* stored backwards: the index finds them */
            hindex[k * 6 + j] = pnt;
            for (int i = 0; i < 3; ++i) hpos[3 * pnt + i] = (float)q[j][i];
        }
    }
    for (int k = 0; k < NP; ++k) {
        poly_corners(q, 5, pc[k], 1.0, ptilt[k], pspin[k], k == 1, 0.0);
        for (int j = 0; j < 5; ++j) for (int i = 0; i < 3; ++i) ppos[3 * (k * 5 + j) + i] = (float)q[j][i];
    }
    {
        double tc[3] = { 6.4, 3.5, -14.5 };
        poly_corners(q, 3, tc, 1.2, 25, 45, false, 0.0);
        for (int j = 0; j < 3; ++j) for (int i = 0; i < 3; ++i) tpos[3 * j + i] = (float)q[j][i];
        const float line[9] = { 5.4f, 0.0f, -15.0f,  7.4f, 0.9f, -14.6f,  6.4f, 0.45f, -14.8f };   /* the third between the others */
        memcpy(tpos + 9, line, sizeof(line));
        double oc[3] = { 7.4f, -3.0f, -16.0f };
        poly_corners(q, 8, oc, 1.1, 30, 75, false, 0.0);
        for (int j = 0; j < 8; ++j) for (int i = 0; i < 3; ++i) opos[3 * j + i] = (float)q[j][i];
    }
    const int hidden_hex = 5;
    uint32_t hvis[2] = { ~0u, ~0u };
    { uint32_t pnt = hindex[hidden_hex * 6 + 2]; hvis[pnt >> 5] &= ~(1u << (pnt & 31)); }
    /* the pentagons' second instance: turned 25 degrees about y, scaled 1.1, behind and aside */
    const double ang = 25 * 3.14159265358979 / 180, sc = 1.1;
    const double tx[3] = { 0.8, 0.5, -3.0 };
    float xf[32] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };
    xf[16] = (float)(sc * cos(ang)); xf[17] = 0; xf[18] = (float)(-sc * sin(ang)); xf[19] = 0;
    xf[20] = 0; xf[21] = (float)sc; xf[22] = 0; xf[23] = 0;
    xf[24] = (float)(sc * sin(ang)); xf[25] = 0; xf[26] = (float)(sc * cos(ang)); xf[27] = 0;
    /* about the pentagons' middle, so the copy stays in view */
    const double mid[3] = { -1.25, -3.35, -15.0 };
    xf[28] = (float)(mid[0] + tx[0] - (xf[16] * mid[0] + xf[24] * mid[2]));
    xf[29] = (float)(mid[1] + tx[1] - xf[21] * mid[1]);
    xf[30] = (float)(mid[2] + tx[2] - (xf[18] * mid[0] + xf[26] * mid[2]));
    xf[31] = 1;

    md_render_batch_t bp[4];
    memset(bp, 0, sizeof(bp));
    bp[0].key = 21; bp[0].kind = MD_RENDER_POLYGONS; bp[0].count = NH; bp[0].points.count = NH * 6;
    bp[0].points.position = f_buf(&fx, hpos, sizeof(hpos));
    bp[0].points.visible[0] = f_buf(&fx, hvis, sizeof(hvis));
    bp[0].index = f_buf(&fx, hindex, sizeof(hindex));
    bp[0].style.polygons.corners = 6; bp[0].style.polygons.radius = 0.22f;
    bp[1].key = 22; bp[1].kind = MD_RENDER_POLYGONS; bp[1].count = NP; bp[1].points.count = NP * 5;
    bp[1].points.position = f_buf(&fx, ppos, sizeof(ppos));
    bp[1].style.polygons.corners = 5; bp[1].style.polygons.radius = 0.15f;
    bp[1].transform = xf; bp[1].transform_count = 2;
    bp[2].key = 23; bp[2].kind = MD_RENDER_POLYGONS; bp[2].count = NT; bp[2].points.count = NT * 3;
    bp[2].points.position = f_buf(&fx, tpos, sizeof(tpos));
    bp[2].style.polygons.corners = 3; bp[2].style.polygons.radius = 0.3f;
    bp[3].key = 24; bp[3].kind = MD_RENDER_POLYGONS; bp[3].count = NO; bp[3].points.count = NO * 8;
    bp[3].points.position = f_buf(&fx, opos, sizeof(opos));
    bp[3].style.polygons.corners = 8; bp[3].style.polygons.radius = 0.2f;
    md_render_frame_t f = frame_of(W, H, cam, bp, 4);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));

    /* the reference, view space (the view is the identity) */
    enum { NREF = NH + 2 * NP + NT + NO };
    poly_ref_t ref[NREF];
    int nref = 0;
    for (int k = 0; k < NH; ++k) {
        if (k == hidden_hex) continue;
        double qq[8][3];
        for (int j = 0; j < 6; ++j) for (int i = 0; i < 3; ++i) qq[j][i] = hpos[3 * hindex[k * 6 + j] + i];
        poly_ref_build(&ref[nref], (const double (*)[3])qq, 6, 0.22f);
        ref[nref].key = 21; ref[nref].linear = (uint32_t)k;
        for (int j = 0; j < 6; ++j) ref[nref].point[j] = hindex[k * 6 + j];
        nref++;
    }
    for (int in = 0; in < 2; ++in) {
        for (int k = 0; k < NP; ++k) {
            double qq[8][3];
            for (int j = 0; j < 5; ++j) {
                const float* p = ppos + 3 * (k * 5 + j);
                const float* m = xf + 16 * in;
                for (int i = 0; i < 3; ++i) qq[j][i] = (double)m[i] * p[0] + (double)m[4 + i] * p[1] + (double)m[8 + i] * p[2] + m[12 + i];
            }
            poly_ref_build(&ref[nref], (const double (*)[3])qq, 5, 0.15f * (in ? sc : 1.0));
            ref[nref].key = 22; ref[nref].linear = (uint32_t)(in * NP + k);
            for (int j = 0; j < 5; ++j) ref[nref].point[j] = (uint32_t)(k * 5 + j);
            nref++;
        }
    }
    for (int k = 0; k < NT; ++k) {
        double qq[8][3];
        for (int j = 0; j < 3; ++j) for (int i = 0; i < 3; ++i) qq[j][i] = tpos[3 * (k * 3 + j) + i];
        poly_ref_build(&ref[nref], (const double (*)[3])qq, 3, 0.3f);
        ref[nref].key = 23; ref[nref].linear = (uint32_t)k;
        for (int j = 0; j < 3; ++j) ref[nref].point[j] = (uint32_t)(k * 3 + j);
        nref++;
    }
    {
        double qq[8][3];
        for (int j = 0; j < 8; ++j) for (int i = 0; i < 3; ++i) qq[j][i] = opos[3 * j + i];
        poly_ref_build(&ref[nref], (const double (*)[3])qq, 8, 0.2f);
        ref[nref].key = 24; ref[nref].linear = 0;
        for (int j = 0; j < 8; ++j) ref[nref].point[j] = (uint32_t)j;
        nref++;
    }
    EXPECT_FALSE(ref[NH - 1 + 2 * NP + 1].plate);         /* the triangle on a line */

    md_render_targets_t tg = md_render_view_targets(fx.view);
    uint32_t* id = malloc(W * H * 8);
    float*    dp = malloc(W * H * 4);
    uint32_t* nb = malloc(W * H * 4);
    ASSERT_TRUE(f_read(&fx, tg.id, id, W * H * 8) && f_read(&fx, tg.depth, dp, W * H * 4) && f_read(&fx, tg.normal, nb, W * H * 4));
    const uint32_t keys[4] = { 21, 22, 23, 24 };
    int bad = 0, covered = 0, judged = 0, nerr = 0;
    double worst = 0;
    double* errs = malloc(W * H * sizeof(double));
    int per_ref[NREF] = { 0 };
    for (uint32_t y = 0; y < H; ++y) {
        for (uint32_t x = 0; x < W; ++x) {
            uint32_t i = y * W + x;
            ray_t r = pixel_ray(&c, W, H, x + 0.5, y + 0.5);
            uint32_t ix = id[2 * i], iy = id[2 * i + 1];
            double clr = polys_ray_clearance(ref, nref, r, 100.0);
            double px15 = 15.0 * 2.0 * c.tan_half / H;
            if (iy == 0) {
                if (clr < -0.5 * px15) { bad++; if (bad < 8) printf("    pixel (%u, %u): empty, a polygon %.3f px inside\n", x, y, -clr / px15); }
                continue;
            }
            covered++;
            if (clr > 0.5 * px15) { bad++; if (bad < 8) printf("    pixel (%u, %u): covered, the polygons %.3f px away\n", x, y, clr / px15); continue; }
            double zv = -c.n / dp[i], tp = zv / r.d[2];
            double P[3] = { r.o[0] + r.d[0] * tp, r.o[1] + r.d[1] * tp, r.o[2] + r.d[2] * tp };
            double pxw = -P[2] * 2.0 * c.tan_half / H;
            int who; double gap;
            double sd = polys_sdf(ref, nref, P, &who, &gap);
            worst = fmax(worst, fabs(sd) / pxw);
            if (fabs(sd) > 0.05 * pxw) { bad++; if (bad < 8) printf("    pixel (%u, %u): %.3f px off the surface\n", x, y, sd / pxw); continue; }
            if (gap > pxw) {
                judged++;
                per_ref[who]++;
                uint32_t slot = (iy >> 20) - 1u;
                if (slot >= 4 || keys[slot] != ref[who].key || ix != ref[who].linear) {
                    bad++;
                    if (bad < 8) printf("    pixel (%u, %u): batch slot %u element %u, expected key %u element %u\n", x, y, slot, ix, ref[who].key, ref[who].linear);
                    continue;
                }
                /* the normal: the distance's gradient */
                double cp[3], nr[3], ng[3];
                poly_ref_dist(&ref[who], P, cp);
                for (int k = 0; k < 3; ++k) nr[k] = P[k] - cp[k];
                double nl = sqrt(dot3(nr, nr));
                if (nl > 0.5 * ref[who].r && fabs(dot3(nr, r.d)) > 0.3 * nl) {
                    for (int k = 0; k < 3; ++k) nr[k] /= nl;
                    oct_decode(nb[i], ng);
                    double cs = dot3(nr, ng);
                    errs[nerr++] = acos(cs > 1 ? 1 : cs < -1 ? -1 : cs) * 180.0 / 3.14159265358979;
                }
            }
            /* nothing in front of it */
            double tt = 0.0;
            for (int it = 0; it < 400 && tt < tp; ++it) {
                double qv[3] = { r.o[0] + r.d[0] * tt, r.o[1] + r.d[1] * tt, r.o[2] + r.d[2] * tt };
                double d = polys_sdf(ref, nref, qv, &who, &gap);
                double qw = -qv[2] * 2.0 * c.tan_half / H;
                if (d < -0.05 * qw) { bad++; if (bad < 8) printf("    pixel (%u, %u): a polygon %.3f px in front\n", x, y, (tp - tt) / qw); break; }
                tt += fmax(fabs(d), 0.05 * qw);
            }
        }
    }
    qsort(errs, (size_t)nerr, sizeof(double), cmp_double);
    double p50 = nerr ? errs[nerr / 2] : 0, p99 = nerr ? errs[(nerr * 99) / 100] : 0;
    printf("    %d pixels covered, %d with a clear polygon; the farthest %.4f px from the solid; normals over %d: median %.3f, 99th %.3f degrees\n",
           covered, judged, worst, nerr, p50, p99);
    EXPECT_EQ(0, bad);
    EXPECT_GT(covered, (int)(W * H / 6));
    EXPECT_GT(judged, covered * 3 / 4);
    for (int k = 0; k < nref; ++k) EXPECT_GT(per_ref[k], 20);  /* every polygon seen, the instance behind and the line too */
    EXPECT_LT(p50, 0.2);
    EXPECT_LT(p99, 1.0);
    free(errs);

    /* picks: at 70% of the way from each polygon's centroid to its second
       corner, what the raster drew there; its attribute the corner nearest
       the hit */
    int picks = 0;
    for (int k = 0; k < nref; ++k) {
        double pq[3];
        for (int i = 0; i < 3; ++i) pq[i] = ref[k].c[i] + 0.7 * (ref[k].q[1][i] - ref[k].c[i]);
        int32_t qx = (int32_t)((pq[0] / -pq[2] / (c.tan_half * c.aspect) * 0.5 + 0.5) * W);
        int32_t qy = (int32_t)((0.5 - pq[1] / -pq[2] / c.tan_half * 0.5) * H);
        if (qx < 0 || qy < 0 || qx >= (int32_t)W || qy >= (int32_t)H) continue;
        uint32_t i = (uint32_t)qy * W + (uint32_t)qx;
        md_render_pick_t pk;
        ASSERT_TRUE(f_pick(&fx, (uint32_t)qx, (uint32_t)qy, &pk));
        if (id[2 * i + 1] == 0) { EXPECT_EQ(0u, (unsigned)pk.key); continue; }
        uint32_t slot = (id[2 * i + 1] >> 20) - 1u;
        ASSERT_LT(slot, 4u);
        EXPECT_EQ(keys[slot], (uint32_t)pk.key);
        EXPECT_EQ(MD_RENDER_POLYGONS, pk.kind);
        EXPECT_EQ(id[2 * i] % bp[slot].count, pk.element);
        EXPECT_EQ(id[2 * i] / bp[slot].count, pk.instance);
        ray_t r = pixel_ray(&c, W, H, qx + 0.5, qy + 0.5);
        double zv = -c.n / dp[i], tp = zv / r.d[2];
        for (int j = 0; j < 3; ++j) EXPECT_NEAR(r.o[j] + r.d[j] * tp, pk.position[j], 1e-4 * 15.0);
        int hit = -1;
        for (int j = 0; j < nref; ++j) if (ref[j].key == keys[slot] && ref[j].linear == id[2 * i]) hit = j;
        ASSERT_GE(hit, 0);
        double P[3] = { pk.position[0], pk.position[1], pk.position[2] }, bd = 1e30;
        uint32_t near = 0;
        for (int j = 0; j < ref[hit].corners; ++j) {
            double e[3] = { P[0] - ref[hit].q[j][0], P[1] - ref[hit].q[j][1], P[2] - ref[hit].q[j][2] };
            if (dot3(e, e) < bd) { bd = dot3(e, e); near = ref[hit].point[j]; }
        }
        EXPECT_EQ(near, pk.attr);
        picks++;
    }
    EXPECT_GT(picks, nref / 2);

    /* state: hexagons 0 and 2 with all their corners, 1 with all but one;
       pentagon 0 (both instances). The mask where the raster shows those. */
    uint32_t hsel[2] = { 0, 0 }, psel = 0x1Fu;
    for (int k = 0; k < 3; ++k) {
        for (int j = 0; j < 6; ++j) {
            if (k == 1 && j == 4) continue;
            uint32_t pnt = hindex[k * 6 + j];
            hsel[pnt >> 5] |= 1u << (pnt & 31);
        }
    }
    bp[0].points.state[0] = f_buf(&fx, hsel, sizeof(hsel));
    bp[1].points.state[0] = f_buf(&fx, &psel, 4);
    f.content[MD_RENDER_CONTENT_STATE]++;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    uint32_t* m = malloc(W * H * 4);
    ASSERT_TRUE(f_read(&fx, tg.mask, m, W * H * 4));
    int mbad = 0, mvis = 0, many = 0;
    for (uint32_t i = 0; i < W * H; ++i) {
        uint32_t slot = id[2 * i + 1] ? (id[2 * i + 1] >> 20) - 1u : 99u, e = id[2 * i];
        bool want = (slot == 0 && (e == 0 || e == 2)) || (slot == 1 && e % NP == 0);
        bool vis = (m[i] & 0xFFu) != 0, any = (m[i] & 0xFF00u) != 0;
        mvis += vis; many += any;
        if (vis != want || (vis && !any) || (m[i] & 0xFFFF0000u)) {
            if (mbad < 6) printf("    pixel (%u, %u): mask %08x, slot %u element %u\n", i % W, i / W, m[i], slot, e);
            mbad++;
        }
    }
    EXPECT_EQ(0, mbad);
    EXPECT_GT(mvis, 300);
    EXPECT_GE(many, mvis);
    printf("    selected: %d pixels in view, %d in all\n", mvis, many);
    free(m);
    bp[0].points.state[0] = 0;
    bp[1].points.state[0] = 0;

    /* shrunk across a fade band: fewer pixels; past it, nothing */
    for (int k = 0; k < 4; ++k) { bp[k].fade.mode = MD_RENDER_FADE_SHRINK; bp[k].fade.out[0] = 14.0f; bp[k].fade.out[1] = 8.0f; }
    f.content[MD_RENDER_CONTENT_GEOMETRY]++;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    ASSERT_TRUE(f_read(&fx, tg.id, id, W * H * 8));
    int shrunk = 0;
    for (uint32_t i = 0; i < W * H; ++i) shrunk += id[2 * i + 1] != 0;
    EXPECT_GT(shrunk, covered / 8);
    EXPECT_LT(shrunk, covered * 3 / 4);
    for (int k = 0; k < 4; ++k) { bp[k].fade.out[0] = 40.0f; bp[k].fade.out[1] = 30.0f; }
    f.content[MD_RENDER_CONTENT_GEOMETRY]++;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    md_gpu_stream_sync(fx.gfx);
    md_render_stats_t st;
    ASSERT_TRUE(md_render_view_stats(fx.view, &st));
    EXPECT_EQ((uint32_t)nref, st.faded);
    ASSERT_TRUE(f_read(&fx, tg.id, id, W * H * 8));
    int left = 0;
    for (uint32_t i = 0; i < W * H; ++i) left += id[2 * i + 1] != 0;
    EXPECT_EQ(0, left);
    printf("    shrunk: %d pixels; past the band: %d\n", shrunk, left);
    for (int k = 0; k < 4; ++k) memset(&bp[k].fade, 0, sizeof(bp[k].fade));

    /* with spheres (the main kernels) and a spline (the same twin) in one
       frame: every pixel as the nearest of the three drawn alone */
    enum { NS = 16, NC = 6 };
    float spos[NS * 3];
    for (int k = 0; k < NS; ++k) {
        const poly_ref_t* p = &ref[(k * 5) % nref];
        for (int i = 0; i < 3; ++i) spos[3 * k + i] = (float)p->q[k % p->corners][i];
        spos[3 * k + 2] += (k % 2) ? 0.5f : -0.5f;
    }
    float cpos[NC * 3];
    uint32_t cidx[NC];
    for (int k = 0; k < NC; ++k) {
        cpos[3 * k + 0] = -7.0f + 2.8f * k; cpos[3 * k + 1] = 1.8f + 0.6f * (float)sin(1.7 * k); cpos[3 * k + 2] = -14.8f + 0.8f * (float)cos(1.1 * k);
        cidx[k] = (uint32_t)k;
    }
    cidx[0] |= MD_RENDER_SPLINE_FIRST; cidx[NC - 1] |= MD_RENDER_SPLINE_LAST;
    md_render_batch_t bs = {0}, bc = {0};
    bs.key = 31; bs.kind = MD_RENDER_SPHERES; bs.count = NS; bs.points.count = NS;
    bs.points.position = f_buf(&fx, spos, sizeof(spos));
    bs.style.spheres.radius_add = 0.4f;
    bs.color.mode = MD_RENDER_COLOR_UNIFORM; bs.color.uniform = 0xFF30C030u;
    bc.key = 32; bc.kind = MD_RENDER_SPLINES; bc.count = NC; bc.points.count = NC;
    bc.points.position = f_buf(&fx, cpos, sizeof(cpos));
    bc.index = f_buf(&fx, cidx, sizeof(cidx));
    bc.style.splines.radius = 0.35f;
    bc.color.mode = MD_RENDER_COLOR_UNIFORM; bc.color.uniform = 0xFFC03030u;
    md_render_batch_t all[6] = { bs, bp[0], bp[1], bp[2], bp[3], bc };
    md_render_view_t vw[4] = { md_render_view_create(fx.r, "spheres"), md_render_view_create(fx.r, "polygons"),
                               md_render_view_create(fx.r, "spline"), md_render_view_create(fx.r, "all") };
    md_render_frame_t fr[4] = { frame_of(W, H, cam, &bs, 1), frame_of(W, H, cam, bp, 4), frame_of(W, H, cam, &bc, 1), frame_of(W, H, cam, all, 6) };
    const uint32_t slot0[3] = { 0, 1, 5 };      /* each one's first slot in `all` */
    for (int k = 0; k < 4; ++k) ASSERT_TRUE(md_render_view_draw(fx.gfx, vw[k], &fr[k], NULL));
    md_gpu_stream_sync(fx.gfx);
    uint32_t* ids[4]; float* dps[4]; uint32_t* cols[4]; uint32_t* nrms[4];
    for (int k = 0; k < 4; ++k) {
        md_render_targets_t t = md_render_view_targets(vw[k]);
        ids[k] = malloc(W * H * 8); dps[k] = malloc(W * H * 4); cols[k] = malloc(W * H * 4); nrms[k] = malloc(W * H * 4);
        ASSERT_TRUE(f_read(&fx, t.id, ids[k], W * H * 8) && f_read(&fx, t.depth, dps[k], W * H * 4) &&
                    f_read(&fx, t.color, cols[k], W * H * 4) && f_read(&fx, t.normal, nrms[k], W * H * 4));
    }
    int mixed_bad = 0, from[3] = { 0, 0, 0 };
    for (uint32_t i = 0; i < W * H; ++i) {
        int w = -1;
        for (int k = 0; k < 3; ++k) if (ids[k][2 * i + 1] != 0 && (w < 0 || dps[k][i] > dps[w][i])) w = k;   /* reverse Z: the larger is nearer */
        if (w < 0) { if (ids[3][2 * i + 1] != 0) mixed_bad++; continue; }
        from[w]++;
        uint32_t ay = ids[w][2 * i + 1];
        uint32_t ey = (ay & 0xFFFFFu) | ((slot0[w] + (ay >> 20)) << 20);
        if (ids[3][2 * i] != ids[w][2 * i] || ids[3][2 * i + 1] != ey || dps[3][i] != dps[w][i] || cols[3][i] != cols[w][i] || nrms[3][i] != nrms[w][i]) {
            if (mixed_bad < 6) printf("    pixel (%u, %u): all %u/%x, alone (%d) %u/%x\n", i % W, i / W, ids[3][2 * i], ids[3][2 * i + 1], w, ids[w][2 * i], ids[w][2 * i + 1]);
            mixed_bad++;
        }
    }
    EXPECT_EQ(0, mixed_bad);
    EXPECT_GT(from[0], 100);
    EXPECT_GT(from[1], 1000);
    EXPECT_GT(from[2], 100);
    printf("    with spheres and a spline: %d, %d and %d pixels\n", from[0], from[1], from[2]);
    for (int k = 0; k < 4; ++k) { free(ids[k]); free(dps[k]); free(cols[k]); free(nrms[k]); md_render_view_destroy(fx.gfx, vw[k]); }
    free(id); free(dp); free(nb);
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

    /* The camera moves a little: a raster that keeps the state lists, and the
       masks still show the selection */
    f.camera.view[12] = 0.05f;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    EXPECT_EQ(MD_RENDER_WORK_RASTER | MD_RENDER_WORK_RESOLVE | MD_RENDER_WORK_MASK, st.work);
    ASSERT_TRUE(f_read(&fx, md_render_view_targets(fx.view).mask, m, W * H * 4));
    #define AT(wx) m[(H / 2) * W + (uint32_t)(((wx) / 4.0 * 0.5 + 0.5) * W)]
    EXPECT_EQ(0u, AT(-1.5));
    EXPECT_EQ(0x0000FF00u | 0x000000FFu, AT(-0.3));
    EXPECT_EQ(0xFF000000u | 0x00FF0000u, AT(2.5));    /* the highlight too */
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

/* ---- Culling and the compute raster ---------------------------------------------
   The options change the cost, never the image: a dense scene of spheres and
   capsules of all sizes (many under a pixel or two, most occluded) through
   the full pipeline must give the ids and depths of the plain hardware raster
   with nothing culled but the frustum. */

typedef struct {
    uint32_t n_sph, n_cap;
    md_gpu_addr_t pos, rad, cap_idx;
    md_render_batch_t b[2];
} dense_t;

static void dense_scene(fix_t* fx, dense_t* d) {
    const int NX = 48, NY = 32, NZ = 12;
    d->n_sph = (uint32_t)(NX * NY * NZ);
    float* pos = malloc(d->n_sph * 12);
    float* rad = malloc(d->n_sph * 4);
    uint32_t seed = 12345;
    #define RND() ((seed = seed * 1664525u + 1013904223u) >> 8) / 16777216.0f
    for (int z = 0, i = 0; z < NZ; ++z) for (int y = 0; y < NY; ++y) for (int x = 0; x < NX; ++x, ++i) {
        pos[3 * i + 0] = (x - NX / 2) * 0.9f + RND() * 0.3f;
        pos[3 * i + 1] = (y - NY / 2) * 0.9f + RND() * 0.3f;
        pos[3 * i + 2] = -8.0f - z * 3.0f - RND();
        rad[i] = z == 0 ? 0.65f : 0.05f + 0.5f * RND() * RND();   /* the front layer is a wall */
    }
    #undef RND
    /* capsules between neighbours along x */
    d->n_cap = 0;
    uint32_t* idx = malloc(d->n_sph * 8);
    for (uint32_t i = 0; i + 1 < d->n_sph; i += 3) { idx[2 * d->n_cap] = i; idx[2 * d->n_cap + 1] = i + 1; d->n_cap++; }
    d->pos = f_buf(fx, pos, d->n_sph * 12);
    d->rad = f_buf(fx, rad, d->n_sph * 4);
    d->cap_idx = f_buf(fx, idx, d->n_cap * 8);
    free(pos); free(rad); free(idx);
    memset(d->b, 0, sizeof(d->b));
    d->b[0].key = 1; d->b[0].kind = MD_RENDER_SPHERES; d->b[0].count = d->n_sph; d->b[0].points.count = d->n_sph;
    d->b[0].points.position = d->pos; d->b[0].points.radius = d->rad;
    d->b[0].style.spheres.radius_scale = 1.0f;
    d->b[1].key = 2; d->b[1].kind = MD_RENDER_CAPSULES; d->b[1].count = d->n_cap; d->b[1].points.count = d->n_sph;
    d->b[1].points.position = d->pos; d->b[1].index = d->cap_idx;
    d->b[1].style.capsules.radius = 0.04f;
}

/* md_render's occlusion settings, for a test of the path its defaults skip
   (in md_render.c, not in its header) */
void md_render_debug_occlusion(float margin_radii, float min_px);

/* ids and depths equal; with depth_only, depths (an index renumbers the elements) */
static int compare_views_ex(fix_t* fx, md_render_view_t a, md_render_view_t b, uint32_t W, uint32_t H, int* covered, bool depth_only) {
    uint32_t* ia = malloc(W * H * 8); uint32_t* ib = malloc(W * H * 8);
    float* da = malloc(W * H * 4); float* db = malloc(W * H * 4);
    int bad = -1;
    if (f_read(fx, md_render_view_targets(a).id, ia, W * H * 8) && f_read(fx, md_render_view_targets(b).id, ib, W * H * 8) &&
        f_read(fx, md_render_view_targets(a).depth, da, W * H * 4) && f_read(fx, md_render_view_targets(b).depth, db, W * H * 4)) {
        bad = 0; *covered = 0;
        for (uint32_t i = 0; i < W * H; ++i) {
            if (ib[2 * i + 1]) (*covered)++;
            if ((!depth_only && (ia[2 * i] != ib[2 * i] || ia[2 * i + 1] != ib[2 * i + 1])) || da[i] != db[i]) {
                if (bad < 8) printf("    pixel (%u, %u): id %u/%x depth %.9g  vs reference %u/%x depth %.9g\n", i % W, i / W,
                                    ia[2 * i], ia[2 * i + 1], da[i], ib[2 * i], ib[2 * i + 1], db[i]);
                bad++;
            }
        }
    }
    free(ia); free(ib); free(da); free(db);
    return bad;
}

static int compare_views(fix_t* fx, md_render_view_t a, md_render_view_t b, uint32_t W, uint32_t H, int* covered) {
    return compare_views_ex(fx, a, b, W, H, covered, false);
}

UTEST(render, culling_and_compute_raster_change_nothing) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    md_gpu_device_info_t info;
    md_gpu_device_info(fx.dev, &info);
    const uint32_t W = 160, H = 120;
    cam_t c = { .tan_half = 0.6f, .aspect = (float)W / H, .n = 0.1f };
    dense_t d;
    dense_scene(&fx, &d);

    md_render_view_t ref = md_render_view_create(fx.r, "reference");
    md_render_frame_t fr = frame_of(W, H, make_camera(&c), d.b, 2);
    fr.options = MD_RENDER_NO_OCCLUSION_CULLING | MD_RENDER_NO_COMPUTE_RASTER;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, ref, &fr, NULL));

    md_render_frame_t f = frame_of(W, H, make_camera(&c), d.b, 2);
    f.small_pixels = 4.0f;
    md_render_stats_t st;
    /* First raster: no pyramid yet, so nothing is occluded in phase 1 */
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    md_gpu_stream_sync(fx.gfx);
    ASSERT_TRUE(md_render_view_stats(fx.view, &st));
    EXPECT_EQ((uint64_t)(d.n_sph + d.n_cap), st.elements);
    EXPECT_EQ(0u, st.occlusion_retested);
    if (info.supports_atomic_int64) EXPECT_GT(st.compute[0], 0u);
    int covered = 0;
    EXPECT_EQ(0, compare_views(&fx, fx.view, ref, W, H, &covered));
    EXPECT_GT(covered, (int)(W * H / 2));

    /* Second raster of the same scene (the positions announced as changed):
       phase 1 culls against the first raster's pyramid */
    f.content[MD_RENDER_CONTENT_GEOMETRY]++;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    md_gpu_stream_sync(fx.gfx);
    ASSERT_TRUE(md_render_view_stats(fx.view, &st));
    EXPECT_EQ(2u, st.raster);
    EXPECT_TRUE(st.refitted);                             /* the geometry was announced as changed */
    EXPECT_GT(st.occluded, (d.n_sph + d.n_cap) / 5);       /* much of a 12-deep stack (elements up to 2 px are drawn, not tested) */
    EXPECT_EQ(st.elements, (uint64_t)st.hidden + st.outside + st.occluded + st.hardware[0] + st.hardware[1] + st.compute[0] + st.compute[1]);
    EXPECT_EQ(0, compare_views(&fx, fx.view, ref, W, H, &covered));
    printf("    elements %llu: outside %u, occluded %u (retested %u), hardware %u + %u, compute %u + %u\n",
           (unsigned long long)st.elements, st.outside, st.occluded, st.occlusion_retested, st.hardware[0], st.hardware[1], st.compute[0], st.compute[1]);

    /* The camera moves: the old pyramid is wrong for it. Phase 2 must recover
       everything phase 1 culled wrongly. */
    md_render_camera_t moved = make_camera(&c);
    moved.view[12] = 2.5f;       /* translate x */
    moved.view[13] = -1.0f;
    f.camera  = moved;
    fr.camera = moved;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, ref, &fr, NULL));
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    md_gpu_stream_sync(fx.gfx);
    ASSERT_TRUE(md_render_view_stats(fx.view, &st));
    EXPECT_GT(st.hardware[1] + st.compute[1], 0u);         /* phase 2 found some */
    EXPECT_FALSE(st.refitted);                             /* only the camera moved: bounds kept */
    EXPECT_EQ(0, compare_views(&fx, fx.view, ref, W, H, &covered));
    printf("    moved:  occluded %u (retested %u), phase 2 drew %u + %u\n", st.occluded, st.occlusion_retested, st.hardware[1], st.compute[1]);

    /* Explicit clusters over the same spheres, compact this time: 4 x 4 tiles
       of each layer, through an index that orders the spheres tile by tile.
       The same image, and the cluster level now culls. */
    {
        const uint32_t NX = 48, NY = 32, NZ = 12;
        uint32_t* perm = malloc(d.n_sph * 4);
        uint32_t* off  = malloc((d.n_sph + 1) * 4);
        uint32_t ne = 0, nc = 0;
        for (uint32_t z = 0; z < NZ; ++z) for (uint32_t ty = 0; ty < NY; ty += 4) for (uint32_t tx = 0; tx < NX; tx += 4) {
            off[nc++] = ne;
            for (uint32_t y = ty; y < ty + 4; ++y) for (uint32_t x = tx; x < tx + 4; ++x) perm[ne++] = (z * NY + y) * NX + x;
        }
        off[nc] = ne;
        md_gpu_addr_t C = f_buf(&fx, off, (nc + 1) * 4);
        md_gpu_addr_t I = f_buf(&fx, perm, ne * 4);
        md_render_batch_t bc[2] = { d.b[0], d.b[1] };
        bc[0].index          = I;
        bc[0].cluster_offset = C;
        bc[0].cluster_count  = nc;
        f.batches = bc;
        ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));     /* new clusters: a refit */
        f.content[MD_RENDER_CONTENT_GEOMETRY]++;
        ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));     /* with the pyramid */
        md_gpu_stream_sync(fx.gfx);
        ASSERT_TRUE(md_render_view_stats(fx.view, &st));
        EXPECT_EQ(nc + (d.n_cap + 63) / 64, st.clusters);
        EXPECT_GT(st.clusters_occluded, nc / 5);           /* a quarter are outside; of the rest, the layers behind the wall */
        EXPECT_EQ(st.elements, (uint64_t)st.hidden + st.outside + st.occluded + st.hardware[0] + st.hardware[1] + st.compute[0] + st.compute[1]);
        EXPECT_EQ(0, compare_views_ex(&fx, fx.view, ref, W, H, &covered, true));
        printf("    tiles:  clusters %u, expanded %u + %u, occluded %u; elements occluded %u\n", st.clusters,
               st.clusters_expanded[0], st.clusters_expanded[1], st.clusters_occluded, st.occluded);
        f.batches = d.b;
        free(off);
        free(perm);
    }

    /* Every element tested for occlusion, the tiny ones too (md_render
       leaves those up to 2 px to its depth test by default): the same image,
       more occluded */
    {
        md_render_debug_occlusion(2.0f, 0.0f);
        md_render_view_t vo = md_render_view_create(fx.r, "every element tested");
        md_render_frame_t fo = frame_of(W, H, make_camera(&c), d.b, 2);
        fo.small_pixels = 4.0f;
        ASSERT_TRUE(md_render_view_draw(fx.gfx, vo, &fo, NULL));
        fo.content[MD_RENDER_CONTENT_GEOMETRY]++;
        ASSERT_TRUE(md_render_view_draw(fx.gfx, vo, &fo, NULL));
        md_gpu_stream_sync(fx.gfx);
        md_render_stats_t so;
        ASSERT_TRUE(md_render_view_stats(vo, &so));
        EXPECT_GT(so.occluded, (d.n_sph + d.n_cap) / 4);
        md_render_view_t ro = md_render_view_create(fx.r, "its reference");
        md_render_frame_t fro = fo;
        fro.options = MD_RENDER_NO_OCCLUSION_CULLING | MD_RENDER_NO_COMPUTE_RASTER;
        ASSERT_TRUE(md_render_view_draw(fx.gfx, ro, &fro, NULL));
        EXPECT_EQ(0, compare_views(&fx, vo, ro, W, H, &covered));
        md_render_view_destroy(fx.gfx, vo);
        md_render_view_destroy(fx.gfx, ro);
        md_render_debug_occlusion(2.0f, 2.0f);
    }

    /* Jittered: the compute raster tests only the pixels the jittered shape
       can reach, so check it against the hardware raster sample by sample */
    {
        md_render_view_t vj = md_render_view_create(fx.r, "jittered");
        md_render_view_t rj = md_render_view_create(fx.r, "jittered reference");
        md_render_frame_t fj = frame_of(W, H, make_camera(&c), d.b, 2);
        fj.jitter = MD_RENDER_JITTER_ALWAYS;
        fj.small_pixels = 4.0f;
        md_render_frame_t fjr = fj;
        fjr.options = MD_RENDER_NO_OCCLUSION_CULLING | MD_RENDER_NO_COMPUTE_RASTER;
        for (int k = 0; k < 5; ++k) {
            ASSERT_TRUE(md_render_view_draw(fx.gfx, vj, &fj, NULL));
            ASSERT_TRUE(md_render_view_draw(fx.gfx, rj, &fjr, NULL));
            EXPECT_EQ(0, compare_views(&fx, vj, rj, W, H, &covered));
        }
        md_render_view_destroy(fx.gfx, vj);
        md_render_view_destroy(fx.gfx, rj);
    }

    /* Selection after culling: the state cull (phase 3) follows a change of
       the bits without a raster, hidden ones included */
    uint32_t words = (d.n_sph + 31) / 32;
    uint32_t* sel = calloc(words, 4);
    const uint32_t deep = d.n_sph - 1;                      /* back of the stack: hidden */
    sel[deep / 32] |= 1u << (deep % 32);
    md_gpu_addr_t S = f_buf(&fx, NULL, words * 4);
    md_gpu_memset(fx.gfx, S, 0, words * 4);
    d.b[0].points.state[0] = S;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));   /* the address is geometry: a raster */
    md_gpu_upload(fx.gfx, S, sel, words * 4);
    f.content[MD_RENDER_CONTENT_STATE]++;
    md_render_status_t rs;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &rs));
    EXPECT_EQ(MD_RENDER_WORK_MASK, rs.work);
    uint32_t* m = malloc(W * H * 4);
    ASSERT_TRUE(f_read(&fx, md_render_view_targets(fx.view).mask, m, W * H * 4));
    int anywhere = 0, visible = 0;
    for (uint32_t i = 0; i < W * H; ++i) { if (m[i] & 0x0000FF00u) anywhere++; if (m[i] & 0xFFu) visible++; }
    EXPECT_GT(anywhere, 0);
    EXPECT_LE(visible, anywhere);
    free(m); free(sel);
    md_render_view_destroy(fx.gfx, ref);
    f_close(&fx);
}

/* Radii from tables: six tables over seven batches (one read by two, one of
   600 entries), types past a table's end among them, flags above the type,
   small spheres for the compute raster and large ones for the hardware. They
   draw what float radii draw, to the bit; a sphere style that ignores radii
   (radius_scale 0, licorice's) draws what no radii do. */
UTEST(render, radius_tables_change_nothing) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    md_gpu_device_info_t info;
    md_gpu_device_info(fx.dev, &info);
    const uint32_t W = 160, H = 120;
    cam_t c = { .tan_half = 0.6f, .aspect = (float)W / H, .n = 0.1f };
    enum { NB = 7, PER = 48, NT = 6 };
    /* the last batch reads table 0 again */
    const uint32_t tsize[NT] = { 8, 600, 13, 5, 300, 7 };
    const uint32_t which[NB] = { 0, 1, 2, 3, 4, 5, 0 };
    float* table[NT];
    md_gpu_addr_t T[NT];
    uint32_t seed = 4711;
    #define RND() ((seed = seed * 1664525u + 1013904223u) >> 8) / 16777216.0f
    for (uint32_t t = 0; t < NT; ++t) {
        table[t] = malloc(tsize[t] * sizeof(float));
        for (uint32_t i = 0; i < tsize[t]; ++i) table[t][i] = 0.05f + 0.6f * RND() * RND();
        T[t] = f_buf(&fx, table[t], tsize[t] * sizeof(float));
    }
    md_render_batch_t bt[NB + 1], bf[NB + 1];         /* by table, by float radius; and one more, licorice's style */
    memset(bt, 0, sizeof(bt));
    memset(bf, 0, sizeof(bf));
    float pts[PER * 4], rad[PER];
    for (uint32_t k = 0; k <= NB; ++k) {
        const uint32_t t = which[k < NB ? k : 0];
        for (uint32_t i = 0; i < PER; ++i) {
            float z = -6.0f - 24.0f * RND();
            pts[4 * i + 0] = (RND() * 2.0f - 1.0f) * 0.55f * -z * c.aspect;
            pts[4 * i + 1] = (RND() * 2.0f - 1.0f) * 0.55f * -z;
            pts[4 * i + 2] = z;
            uint32_t type = i % 11 == 10 ? tsize[t] + i : (uint32_t)(RND() * tsize[t]) % tsize[t];   /* some past the end: radius 0 */
            uint32_t word = type | ((i * 2654435761u) & 0xFFFF0000u);                           /* flags above, ignored */
            memcpy(&pts[4 * i + 3], &word, 4);
            rad[i] = type < tsize[t] ? table[t][type] : 0.0f;
        }
        md_gpu_addr_t P = f_buf(&fx, pts, sizeof(pts));
        md_render_batch_t b = {0};
        b.key = 10 + k; b.kind = MD_RENDER_SPHERES; b.count = PER; b.points.count = PER;
        b.points.position = P; b.points.stride = 16;
        b.style.spheres.radius_scale = 1.0f;
        bt[k] = bf[k] = b;
        bt[k].points.radius = P + 12; bt[k].points.radius_stride = 16;
        bt[k].points.radius_table = T[t]; bt[k].points.radius_table_count = tsize[t];
        if (k < NB) {
            bf[k].points.radius = f_buf(&fx, rad, sizeof(rad));
        } else {
            bt[k].style.spheres.radius_scale = bf[k].style.spheres.radius_scale = 0.0f;
            bt[k].style.spheres.radius_add   = bf[k].style.spheres.radius_add   = 0.12f;
        }
    }
    #undef RND

    md_render_view_t ref = md_render_view_create(fx.r, "float radii");
    md_render_frame_t fr = frame_of(W, H, make_camera(&c), bf, NB + 1);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, ref, &fr, NULL));

    md_render_frame_t f = frame_of(W, H, make_camera(&c), bt, NB + 1);
    md_render_stats_t st;
    int covered = 0;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    md_gpu_stream_sync(fx.gfx);
    ASSERT_TRUE(md_render_view_stats(fx.view, &st));
    EXPECT_GT(st.hardware[0] + st.hardware[1], 0u);
    if (info.supports_atomic_int64) EXPECT_GT(st.compute[0] + st.compute[1], 0u);
    EXPECT_EQ(0, compare_views(&fx, fx.view, ref, W, H, &covered));
    EXPECT_GT(covered, (int)(W * H / 20));
    /* again, the positions announced as changed: a refit, and phase 1
       against the first raster's pyramid */
    f.content[MD_RENDER_CONTENT_GEOMETRY]++;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    EXPECT_EQ(0, compare_views(&fx, fx.view, ref, W, H, &covered));
    printf("    %u spheres in %u batches, %u tables: %d pixels covered, hardware %u + %u, compute %u + %u\n",
           (NB + 1) * PER, NB + 1, NT, covered, st.hardware[0], st.hardware[1], st.compute[0], st.compute[1]);

    md_render_view_destroy(fx.gfx, ref);
    for (uint32_t t = 0; t < NT; ++t) free(table[t]);
    f_close(&fx);
}

/* ---- Fading with distance ------------------------------------------------- */

/* md_render.c's fade_band and md_render.slang's fade_t and fade_keeps,
   replicated: a band {s0, s1} for this camera as t = saturate(a d + b) */
static void fade_ab(const float band[2], float px, bool ortho, float ab[2]) {
    const double d0 = 1.0 / band[0], d1 = 1.0 / band[1];
    if (ortho) { ab[0] = 0.0f; ab[1] = (float)((1.0 / px - d0) / (d1 - d0)); return; }
    const double a = 1.0 / (px * (d1 - d0));
    ab[0] = (float)a;
    ab[1] = (float)(-px * d0 * a);
}
static float fade_t_cpu(const float ab[2], float d) {
    float t = fmaf(ab[0], d, ab[1]);
    return t < 0.0f ? 0.0f : t > 1.0f ? 1.0f : t;
}
static float fade_h_cpu(uint32_t x) {
    x ^= x >> 16; x *= 0x7FEB352Du;
    x ^= x >> 15; x *= 0x846CA68Bu;
    x ^= x >> 16;
    return (float)(x >> 8) * (1.0f / 16777216.0f);
}

typedef struct fade_scene_t {
    uint32_t      n;
    float*        z;               /* view-space z of each sphere, camera at rest */
    md_gpu_addr_t vis[2];          /* the reference's visibility bits, batch 0 and 1 */
    md_render_batch_t b[2];        /* fading */
    md_render_batch_t r[2];        /* the reference: no fade, visibility bits */
} fade_scene_t;

/* The reference's bits for a THIN handover across `band` (batch 0 fading
   out, 1 in), the camera moved by dz along view z; how near the nearest
   sphere's hash is to its t */
static double fade_split(fix_t* fx, fade_scene_t* fs, const float band[2], float px, float dz, bool ortho) {
    uint32_t words = (fs->n + 31) / 32;
    uint32_t* bits = calloc(2 * words, 4);
    float ab[2];
    fade_ab(band, px, ortho, ab);
    double margin = 1e30;
    for (uint32_t i = 0; i < fs->n; ++i) {
        float d = -(fs->z[i] + dz);
        float t = fade_t_cpu(ab, d), h = fade_h_cpu(i);
        uint32_t k = h >= t ? 0u : 1u;
        bits[k * words + i / 32] |= 1u << (i % 32);
        if (t > 0.0f && t < 1.0f && fabs((double)h - t) < margin) margin = fabs((double)h - t);
    }
    md_gpu_upload(fx->gfx, fs->vis[0], bits, words * 4);
    md_gpu_upload(fx->gfx, fs->vis[1], bits + words, words * 4);
    free(bits);
    return margin;
}

static int compare_masks(fix_t* fx, md_render_view_t a, md_render_view_t b, uint32_t W, uint32_t H, int* set) {
    uint32_t* ma = malloc(W * H * 4); uint32_t* mb = malloc(W * H * 4);
    int bad = -1;
    if (f_read(fx, md_render_view_targets(a).mask, ma, W * H * 4) && f_read(fx, md_render_view_targets(b).mask, mb, W * H * 4)) {
        bad = 0; *set = 0;
        for (uint32_t i = 0; i < W * H; ++i) { if (mb[i]) (*set)++; if (ma[i] != mb[i]) bad++; }
    }
    free(ma); free(mb);
    return bad;
}

/* THIN: one batch thinning out across a band, another thinning in across the
   same band, over the same spheres, draw what the two draw with each sphere
   given by visibility bits to the side its hash puts it on: ids, depths and
   masks to the bit, and each sphere counted faded exactly once. The spheres
   are in depth order, so clusters of 64 leave by their bounds and those in
   the band are decided sphere by sphere. A camera move changes who draws
   what without a refit; orthographic, every sphere is equally far through. */
UTEST(render, fade_thin_hands_over) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 160, H = 120;
    cam_t c = { .tan_half = 0.6f, .aspect = (float)W / H, .n = 0.1f };
    enum { N = 512 };
    fade_scene_t fs;
    memset(&fs, 0, sizeof(fs));
    fs.n = N;
    fs.z = malloc(N * sizeof(float));
    float* pos = malloc(N * 12);
    float* rad = malloc(N * 4);
    uint32_t sel[N / 32];
    for (uint32_t i = 0; i < N / 32; ++i) sel[i] = 0x55555555u;   /* every other sphere selected: masks */
    uint32_t seed = 99;
    #define RND() ((seed = seed * 1664525u + 1013904223u) >> 8) / 16777216.0f
    for (uint32_t i = 0; i < N; ++i) {          /* farther with i: clusters of 64 by depth */
        float z = -4.0f - 56.0f * (float)i / N;
        pos[3 * i + 0] = (RND() * 2.0f - 1.0f) * 0.5f * -z * c.tan_half * c.aspect;
        pos[3 * i + 1] = (RND() * 2.0f - 1.0f) * 0.5f * -z * c.tan_half;
        pos[3 * i + 2] = z;
        rad[i] = 0.05f + 0.15f * RND();
        fs.z[i] = z;
    }
    #undef RND
    md_gpu_addr_t P = f_buf(&fx, pos, N * 12), R = f_buf(&fx, rad, N * 4), S = f_buf(&fx, sel, sizeof(sel));
    fs.vis[0] = f_buf(&fx, NULL, N / 8);
    fs.vis[1] = f_buf(&fx, NULL, N / 8);
    md_render_camera_t cam = make_camera(&c);
    const float px = fabsf(cam.proj[5]) * 0.5f * (float)H;
    /* the band: from sphere 150's depth to sphere 300's, inside clusters */
    float band[2] = { px / -fs.z[150], px / -fs.z[300] };
    for (int k = 0; k < 2; ++k) {
        md_render_batch_t b = {0};
        b.key = 1 + k; b.kind = MD_RENDER_SPHERES; b.count = N; b.points.count = N;
        b.points.position = P; b.points.radius = R; b.points.state[0] = S;
        b.style.spheres.radius_scale = k ? 1.5f : 1.0f;         /* the stand-in: larger */
        b.color.mode = MD_RENDER_COLOR_UNIFORM;
        b.color.uniform = k ? 0xFFFF0000u : 0xFF0000FFu;
        fs.b[k] = fs.r[k] = b;
        fs.r[k].points.visible[0] = fs.vis[k];
        fs.b[k].fade.mode = MD_RENDER_FADE_THIN;
        memcpy(k ? fs.b[k].fade.in : fs.b[k].fade.out, band, sizeof(band));
    }
    while (fade_split(&fx, &fs, band, px, 0.0f, false) < 1e-4) band[1] *= 0.9993f;
    memcpy(fs.b[0].fade.out, band, sizeof(band));
    memcpy(fs.b[1].fade.in,  band, sizeof(band));

    md_render_view_t ref = md_render_view_create(fx.r, "visibility bits");
    md_render_frame_t f  = frame_of(W, H, cam, fs.b, 2);
    md_render_frame_t fr = frame_of(W, H, cam, fs.r, 2);
    md_render_stats_t st;
    int covered = 0, set = 0;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, ref, &fr, NULL));
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    md_gpu_stream_sync(fx.gfx);
    ASSERT_TRUE(md_render_view_stats(fx.view, &st));
    EXPECT_EQ((uint32_t)N, st.faded);
    EXPECT_EQ(st.elements, (uint64_t)st.hidden + st.outside + st.faded + st.occluded + st.hardware[0] + st.hardware[1] + st.compute[0] + st.compute[1]);
    EXPECT_LT(st.clusters_expanded[0] + st.clusters_expanded[1], 16u);     /* 8 clusters a batch: some left whole */
    EXPECT_EQ(0, compare_views(&fx, fx.view, ref, W, H, &covered));
    EXPECT_GT(covered, (int)(W * H / 20));
    EXPECT_EQ(0, compare_masks(&fx, fx.view, ref, W, H, &set));
    EXPECT_GT(set, 0);
    printf("    %u spheres, band %.2f .. %.2f px per unit: %d pixels covered, %d masked; faded %u, clusters expanded %u\n",
           (unsigned)N, band[0], band[1], covered, set, st.faded, st.clusters_expanded[0] + st.clusters_expanded[1]);

    /* The camera backs off: more spheres to the stand-in, the bounds kept */
    float dz = 6.0f;
    while (fade_split(&fx, &fs, band, px, dz, false) < 1e-4) dz += 0.0137f;
    f.camera.view[14] = fr.camera.view[14] = dz;
    fr.content[MD_RENDER_CONTENT_GEOMETRY]++;           /* the reference's bits changed */
    ASSERT_TRUE(md_render_view_draw(fx.gfx, ref, &fr, NULL));
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    md_gpu_stream_sync(fx.gfx);
    ASSERT_TRUE(md_render_view_stats(fx.view, &st));
    EXPECT_FALSE(st.refitted);
    EXPECT_EQ((uint32_t)N, st.faded);
    EXPECT_EQ(0, compare_views(&fx, fx.view, ref, W, H, &covered));
    EXPECT_EQ(0, compare_masks(&fx, fx.view, ref, W, H, &set));

    /* Orthographic: every sphere equally far through the band */
    cam_t co = { .ortho = true, .half_h = 6.0f, .aspect = (float)W / H, .n = 0.0f, .far_ = 80.0f };
    md_render_camera_t camo = make_camera(&co);
    const float pxo = fabsf(camo.proj[5]) * 0.5f * (float)H;
    band[0] = pxo * 1.5f; band[1] = pxo * 0.5f;                /* halfway through, near enough */
    while (fade_split(&fx, &fs, band, pxo, 0.0f, true) < 1e-4) band[1] *= 0.9993f;
    memcpy(fs.b[0].fade.out, band, sizeof(band));
    memcpy(fs.b[1].fade.in,  band, sizeof(band));
    f.camera = fr.camera = camo;
    fr.content[MD_RENDER_CONTENT_GEOMETRY]++;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, ref, &fr, NULL));
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    md_gpu_stream_sync(fx.gfx);
    ASSERT_TRUE(md_render_view_stats(fx.view, &st));
    EXPECT_EQ((uint32_t)N, st.faded);
    EXPECT_EQ(0, compare_views(&fx, fx.view, ref, W, H, &covered));
    EXPECT_GT(covered, 0);

    md_render_view_destroy(fx.gfx, ref);
    free(pos); free(rad); free(fs.z);
    f_close(&fx);
}

/* SHRINK: a sphere's radius scales down to 0 across `out` and up from 0
   across `in`, linearly in depth. Against analytic spheres of those radii
   (coverage, depth, a pick on one in the band), and the spheres past a band
   counted faded. */
UTEST(render, fade_shrink_scales_radii) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 128, H = 96;
    cam_t c = { .tan_half = 0.5f, .aspect = (float)W / H, .n = 0.1f };
    md_render_camera_t cam = make_camera(&c);
    const float px = fabsf(cam.proj[5]) * 0.5f * (float)H;
    enum { N = 12 };
    double ctr[N][3], rexp[N];
    float pos[N * 3], rad[N];
    /* in: nearer than depth 6 absent, full from 8; out: full to 12, gone at 16 */
    md_render_batch_t b = {0};
    b.key = 5; b.kind = MD_RENDER_SPHERES; b.count = N; b.points.count = N;
    b.style.spheres.radius_scale = 1.0f;
    b.fade.mode = MD_RENDER_FADE_SHRINK;
    b.fade.in[0]  = px / 6.0f;  b.fade.in[1]  = px / 8.0f;
    b.fade.out[0] = px / 12.0f; b.fade.out[1] = px / 16.0f;
    float ai[2], ao[2];
    fade_ab(b.fade.in, px, false, ai);
    fade_ab(b.fade.out, px, false, ao);
    uint32_t gone = 0;
    for (int i = 0; i < N; ++i) {
        double z = -(5.0 + 12.0 * i / (N - 1));                     /* 5 .. 17 */
        ctr[i][0] = (i % 4 - 1.5) * 0.25 * -z * c.tan_half * c.aspect;
        ctr[i][1] = ((i / 4) % 3 - 1.0) * 0.3 * -z * c.tan_half;
        ctr[i][2] = z;
        rad[i] = 0.12f * (float)-z * c.tan_half;                    /* about the same on screen */
        for (int k = 0; k < 3; ++k) pos[3 * i + k] = (float)ctr[i][k];
        float d = (float)-z;
        double f = (1.0 - fade_t_cpu(ao, d)) * fade_t_cpu(ai, d);
        rexp[i] = rad[i] * f;
        if (rexp[i] <= 0.0) gone++;
    }
    b.points.position = f_buf(&fx, pos, sizeof(pos));
    b.points.radius   = f_buf(&fx, rad, sizeof(rad));
    md_render_frame_t f = frame_of(W, H, cam, &b, 1);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    md_gpu_stream_sync(fx.gfx);
    int judged = 0;
    EXPECT_EQ(0, check_spheres(&fx, &c, cam, W, H, (const double (*)[3])ctr, rexp, N, &judged));
    EXPECT_GT(judged, (int)(W * H * 9 / 10));
    md_render_stats_t st;
    ASSERT_TRUE(md_render_view_stats(fx.view, &st));
    EXPECT_GT(gone, 0u);
    EXPECT_EQ(gone, st.faded);
    /* a pick on a sphere in the out band: the hit on the shrunk sphere */
    int pi = -1;
    for (int i = 0; i < N; ++i) if (rexp[i] > 0.0 && rexp[i] < rad[i] * 0.9) { pi = i; if (ctr[i][2] < -12.5) break; }
    ASSERT_GE(pi, 0);
    uint32_t qx = (uint32_t)((ctr[pi][0] / -ctr[pi][2] / (c.tan_half * c.aspect) * 0.5 + 0.5) * W);
    uint32_t qy = (uint32_t)((0.5 - ctr[pi][1] / -ctr[pi][2] / c.tan_half * 0.5) * H);
    md_render_pick_t pk;
    ASSERT_TRUE(f_pick(&fx, qx, qy, &pk));
    EXPECT_EQ((uint32_t)pi, pk.element);
    ray_t r = pixel_ray(&c, W, H, qx + 0.5, qy + 0.5);
    double dist, tt = sphere_t(r, ctr[pi], rexp[pi], &dist);
    for (int k = 0; k < 3; ++k) EXPECT_NEAR(r.o[k] + r.d[k] * tt, pk.position[k], 1e-3);
    printf("    %d spheres from depth 5 to 17, in 6 -> 8, out 12 -> 16: %u gone, %d pixels judged; picked %d (radius %.3f of %.3f)\n",
           N, gone, judged, pi, rexp[pi], rad[pi]);
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
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, md_str_from_cstr(ethanol_pdb), NULL));
    ASSERT_TRUE(md_system_infer(&sys, &state, MD_SYSTEM_INFER_ALL, NULL));
    ASSERT_EQ(8u, (uint32_t)(sys.bond.count + sys.bond.noncovalent_count));

    uint64_t content[MD_RENDER_CONTENT_COUNT] = {0};

    /* Topology: plan, the caller allocates, write */
    md_temp_scope_t temp = md_temp_begin();
    md_render_mol_plan_t plan;
    ASSERT_TRUE(md_render_mol_plan(&plan, &sys, NULL, temp));
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
    md_render_mol_t mol = md_render_mol_create(fx.dev, NULL);
    ASSERT_TRUE(mol != NULL);
    md_gpu_addr_t dyn = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, topo.dynamic_size).gpu;
    ASSERT_TRUE(md_render_mol_update(mol, fx.gfx, &topo, dyn, &atoms, 0, content));   /* the copy that is drawn */

    md_render_mol_rep_t rep = {0};
    rep.key  = 42;
    rep.type = MD_RENDER_MOL_REP_LICORICE;
    rep.params.licorice.radius = 0.15f;
    rep.atom_color = f_buf(&fx, colors, sizeof(colors));
    rep.color.mode = MD_RENDER_COLOR_POINTS;

    md_render_batch_t batches[8];
    uint32_t nb = md_render_mol_emit(&topo, dyn, &atoms, &rep, batches, 8);
    ASSERT_EQ(2u, nb);                          /* atoms, covalent bonds */
    EXPECT_EQ(MD_RENDER_CAPSULES, batches[0].kind);
    EXPECT_EQ(MD_RENDER_SPHERES, batches[1].kind);
    EXPECT_EQ(8u, batches[0].count);
    /* emit is pure: asking with no room says how many it needs */
    EXPECT_EQ(2u, md_render_mol_emit(&topo, dyn, &atoms, &rep, NULL, 0));
    /* a representation's fade goes on every batch it emits */
    {
        md_render_mol_rep_t rw = rep;
        rw.fade.mode = MD_RENDER_FADE_THIN;
        rw.fade.out[0] = 2.5f; rw.fade.out[1] = 1.25f;
        md_render_batch_t bw[8];
        ASSERT_EQ(2u, md_render_mol_emit(&topo, dyn, &atoms, &rw, bw, 8));
        for (int k = 0; k < 2; ++k) { EXPECT_EQ((uint32_t)MD_RENDER_FADE_THIN, bw[k].fade.mode); EXPECT_EQ(2.5f, bw[k].fade.out[0]); EXPECT_EQ(1.25f, bw[k].fade.out[1]); }
        EXPECT_EQ(0u, batches[0].fade.mode);
    }

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
    md_gpu_free(fx.gfx, dyn);
    md_gpu_free(fx.gfx, topo.base);
    md_render_mol_destroy(mol);
    md_system_free(&sys);
    f_close(&fx);
}

/* ---- Clusters of the molecular layer -------------------------------------------
   The topology block holds every atom and bond once, in cluster order: static
   clusters of whole residues (at most 64 atoms, one molecule each), then one
   per grid cell holding whole small molecules. Checked on the host image of a
   block (the plan's, or a copy of the GPU's). */

static int32_t* molecule_of_atoms(const md_system_t* sys) {
    int32_t* mol = malloc((sys->atom.count + 1) * sizeof(int32_t));
    for (size_t i = 0; i < sys->atom.count; ++i) mol[i] = -1;
    for (size_t m = 0; m < sys->molecule.count; ++m)
        for (uint32_t k = sys->molecule.offset[m]; k < sys->molecule.offset[m + 1]; ++k) mol[sys->molecule.atom_idx[k]] = (int32_t)m;
    return mol;
}

/* A topology image against its system and resident set (NULL: all): every
   resident atom at one point, no other atom at any; static clusters whole
   residues of one molecule; a small molecule's resident atoms contiguous;
   exactly the bonds with both ends resident, between their atoms' points */
static int check_topo_image(const md_system_t* sys, const md_render_mol_topo_t* T, const uint8_t* img, const md_bitset_t* resident) {
    int bad = 0;
    #define FAIL(...) do { if (bad++ < 6) { printf("    "); printf(__VA_ARGS__); printf("\n"); } } while (0)
    const uint32_t n = T->atom_count, np = T->point_count, nb = T->bond_count;
    const uint32_t nb_sys = (uint32_t)(sys->bond.count + sys->bond.noncovalent_count);
    #define RES(a) (!resident || md_bitset_test(*resident, (a)))
    const uint32_t* order = (const uint32_t*)(img + T->atom_order);
    const uint32_t* coff  = (const uint32_t*)(img + T->atom_cluster_offset);
    int32_t* mol = molecule_of_atoms(sys);
    uint8_t* seen = calloc(n + nb_sys + 1, 1);
    uint32_t n_res = 0;
    for (uint32_t a = 0; a < n; ++a) n_res += RES(a);
    if (np != n_res) FAIL("%u points for %u resident atoms", np, n_res);
    const uint32_t n_static = T->atom_cluster_count - T->grid_cells;
    if (coff[0] != 0 || coff[T->atom_cluster_count] != np) FAIL("atom cluster offsets span %u..%u of %u", coff[0], coff[T->atom_cluster_count], np);
    for (uint32_t c = 0; c < T->atom_cluster_count; ++c) {
        if (coff[c + 1] < coff[c]) { FAIL("atom cluster %u decreasing", c); continue; }
        uint32_t size = coff[c + 1] - coff[c];
        bool dynamic = c >= n_static;
        if (!dynamic && (size == 0 || size > 64)) FAIL("static cluster %u has %u atoms", c, size);
        for (uint32_t k = coff[c]; k < coff[c + 1] && k < np; ++k) {
            uint32_t a = order[k];
            if (a >= n || seen[a]) { FAIL("atom %u at %u: out of range or twice", a, k); continue; }
            if (!RES(a)) FAIL("atom %u at %u is not resident", a, k);
            seen[a] = 1;
            bool small = mol[a] >= 0 && sys->molecule.offset[mol[a] + 1] - sys->molecule.offset[mol[a]] <= MD_RENDER_MOL_SMALL;
            if (small != dynamic) FAIL("atom %u (small %d) in a %s cluster", a, small, dynamic ? "grid" : "static");
            if (!dynamic && mol[a] != mol[order[coff[c]]]) FAIL("static cluster %u crosses molecules", c);
            /* a small molecule is contiguous: the slots after its first resident atom hold the rest */
            if (dynamic && mol[a] >= 0) {
                uint32_t mb = sys->molecule.offset[mol[a]], me = sys->molecule.offset[mol[a] + 1];
                uint32_t rest[MD_RENDER_MOL_SMALL], nr = 0;
                for (uint32_t j = mb; j < me; ++j) if (RES((uint32_t)sys->molecule.atom_idx[j])) rest[nr++] = (uint32_t)sys->molecule.atom_idx[j];
                if (nr && rest[0] == a) {
                    if (k + nr > coff[c + 1]) FAIL("molecule %d split by cluster %u", mol[a], c);
                    for (uint32_t j = 1; j < nr; ++j) if (order[k + j] != rest[j]) { FAIL("molecule %d not contiguous", mol[a]); break; }
                }
            }
        }
    }
    /* the internal map from atom to point: the inverse for the static atoms
       (a small molecule's are kept only with mixed bonds; those show in the
       bonds' points below) */
    const uint32_t* point = (const uint32_t*)(img + T->mixed_point);
    const uint32_t n_static_pts = np - T->small_atom_count;
    for (uint32_t p = 0; p < n_static_pts; ++p) if (order[p] < n && point[order[p]] != p) { FAIL("mixed_point is not the inverse at static point %u", p); break; }
    for (uint32_t a = 0; a < n; ++a) if (!RES(a) && point[a] != UINT32_MAX) { FAIL("atom %u is not resident but has point %u", a, point[a]); break; }
    /* bonds are between points: atom order[p] is at point p */
    const uint32_t* pairs = (const uint32_t*)(img + T->bond_pairs);
    const uint32_t* ids   = (const uint32_t*)(img + T->bond_ids);
    memset(seen, 0, n + nb_sys + 1);
    for (uint32_t j = 0; j < nb; ++j) {
        uint32_t i = ids[j];
        if (i >= nb_sys || seen[i]) { FAIL("bond id %u at %u: out of range or twice", i, j); continue; }
        seen[i] = 1;
        if (pairs[2 * j] >= np || pairs[2 * j + 1] >= np ||
            order[pairs[2 * j]] != (uint32_t)sys->bond.pairs[i].idx[0] || order[pairs[2 * j + 1]] != (uint32_t)sys->bond.pairs[i].idx[1]) FAIL("bond %u: points do not match its id", j);
    }
    for (uint32_t i = 0; i < nb_sys; ++i) {
        uint32_t a = (uint32_t)sys->bond.pairs[i].idx[0], b = (uint32_t)sys->bond.pairs[i].idx[1];
        bool expect = a < n && b < n && RES(a) && RES(b);
        if (expect != (bool)seen[i]) { FAIL("bond %u (%u-%u): %s", i, a, b, expect ? "missing" : "kept with an end not resident"); }
    }
    for (uint32_t c = 0; c < MD_RENDER_MOL_BOND_CLASS_COUNT; ++c) {
        const uint32_t* bo = (const uint32_t*)(img + T->bond_cluster_offset[c]);
        uint32_t size = T->bond_class_offset[c + 1] - T->bond_class_offset[c], cc = T->bond_cluster_count[c];
        if (bo[0] != 0 || bo[cc] != size) FAIL("class %u: clusters span %u..%u of %u", c, bo[0], bo[cc], size);
        uint32_t n_grid = (c == 0 && T->small_bond_count) ? T->grid_cells : 0;
        for (uint32_t k = 0; k < cc; ++k) {
            if (bo[k + 1] < bo[k]) FAIL("class %u cluster %u decreasing", c, k);
            else if (k < cc - n_grid && bo[k + 1] - bo[k] > 64) FAIL("class %u cluster %u has %u bonds", c, k, bo[k + 1] - bo[k]);
        }
    }
    #undef RES
    #undef FAIL
    free(mol); free(seen);
    return bad;
}

static bool load_system(md_system_t* sys, md_system_state_t* state, const char* file) {
    md_str_t path = md_str_from_cstr(file);
    bool ok = md_str_ends_with(path, MD_STR_LIT(".gro")) ? md_gro_system_init_from_file(sys, state, path, NULL)
            : md_pdb_system_init_from_file(sys, state, path, NULL);
    return ok && md_system_infer(sys, state, MD_SYSTEM_INFER_ALL, NULL);
}

static uint8_t* read_block(fix_t* fx, md_gpu_addr_t base, uint64_t size);

/* The plan of systems with chains, ligands and small molecules; then (with a
   GPU) the same after a regroup, which moves the small molecules and with
   them the bonds that reach them from elsewhere (mixed) */
UTEST(render, mol_plan_clusters) {
    const char* files[] = { MD_UNITTEST_DATA_DIR "/tubulin-A-B.pdb", MD_UNITTEST_DATA_DIR "/1a64.pdb", MD_UNITTEST_DATA_DIR "/tryptophan-md.gro" };
    fix_t fx;
    bool gpu = f_open(&fx);
    md_render_mol_t mol = gpu ? md_render_mol_create(fx.dev, NULL) : NULL;
    for (size_t f = 0; f < sizeof(files) / sizeof(files[0]); ++f) {
        md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
        md_system_t sys = { .alloc = alloc };
        md_system_state_t state = { .alloc = alloc };
        ASSERT_TRUE(load_system(&sys, &state, files[f]));
        md_temp_scope_t temp = md_temp_begin();
        md_render_mol_plan_t plan;
        ASSERT_TRUE(md_render_mol_plan(&plan, &sys, NULL, temp));
        const md_render_mol_topo_t* T = &plan.layout;
        EXPECT_EQ(0, check_topo_image(&sys, T, plan.data, NULL));
        uint32_t n_static = T->atom_count - T->small_atom_count, cl_static = T->atom_cluster_count - T->grid_cells;
        printf("    %s: %u atoms: %u in %u residue clusters (%.1f per cluster), %u in %u small molecules (%u kinds) over %u cells; %u bonds (%u mixed)\n",
               strrchr(files[f], '/') + 1, T->atom_count, n_static, cl_static, cl_static ? (double)n_static / cl_static : 0.0,
               T->small_atom_count, T->small_molecule_count, T->kind_count, T->grid_cells, T->bond_count, T->mixed_bond_count);
        EXPECT_LE(T->kind_count, T->small_molecule_count);
        if (n_static > 1000) EXPECT_GT((double)n_static / cl_static, 40.0);   /* residues pack well */
        if (mol) {
            md_render_mol_topo_t topo;
            ASSERT_TRUE(md_render_mol_write(fx.gfx, &topo, &plan, md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, plan.size).gpu, NULL));
            md_gpu_addr_t P   = f_buf(&fx, state.xyz, T->atom_count * 12);
            md_gpu_addr_t dyn = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, topo.dynamic_size).gpu;
            md_render_mol_atoms_t atoms = {0};
            atoms.position = P;
            md_render_mol_topo_t rel = topo;
            rel.atom_order -= topo.base; rel.mixed_point -= topo.base; rel.atom_cluster_offset -= topo.base;
            rel.bond_pairs -= topo.base; rel.bond_ids -= topo.base;
            for (int c = 0; c < MD_RENDER_MOL_BOND_CLASS_COUNT; ++c) rel.bond_cluster_offset[c] -= topo.base;
            /* twice: the second from the first's records (the other buffer),
               mixed bonds and all */
            for (int rg = 0; rg < 2; ++rg) {
                ASSERT_TRUE(md_render_mol_update(mol, fx.gfx, &topo, dyn, &atoms, MD_RENDER_MOL_UPDATE_REGROUP, NULL));
                uint8_t* img = read_block(&fx, topo.base, topo.size);
                EXPECT_EQ(0, check_topo_image(&sys, &rel, img, NULL));
                free(img);
            }
            md_gpu_free(fx.gfx, P); md_gpu_free(fx.gfx, dyn); md_gpu_free(fx.gfx, topo.base);
        }
        md_temp_end(temp);
        md_system_free(&sys);
        md_vm_arena_destroy(alloc);
    }
    if (mol) md_render_mol_destroy(mol);
    if (gpu) f_close(&fx);
}

static uint8_t* read_block(fix_t* fx, md_gpu_addr_t base, uint64_t size) {
    md_gpu_mem_t rb = md_gpu_malloc(fx->gfx, MD_GPU_MEM_HOST_READ, size);
    if (!rb.cpu) return NULL;
    md_gpu_copy(fx->gfx, rb.gpu, base, size);
    md_gpu_stream_sync(fx->gfx);
    uint8_t* out = malloc(size);
    memcpy(out, rb.cpu, size);
    md_gpu_free(fx->gfx, rb.gpu);
    return out;
}

/* Mean bounding-box diagonal of the grid clusters */
static double grid_cluster_extent(const md_render_mol_topo_t* T, const uint8_t* img, const float* xyz) {
    const uint32_t* order = (const uint32_t*)(img + T->atom_order);
    const uint32_t* coff  = (const uint32_t*)(img + T->atom_cluster_offset);
    double sum = 0; uint32_t cnt = 0;
    for (uint32_t c = T->atom_cluster_count - T->grid_cells; c < T->atom_cluster_count; ++c) {
        if (coff[c + 1] == coff[c]) continue;
        double lo[3] = { 1e30, 1e30, 1e30 }, hi[3] = { -1e30, -1e30, -1e30 };
        for (uint32_t k = coff[c]; k < coff[c + 1]; ++k)
            for (int j = 0; j < 3; ++j) { double v = xyz[3 * order[k] + j]; if (v < lo[j]) lo[j] = v; if (v > hi[j]) hi[j] = v; }
        sum += sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) + (hi[2] - lo[2]) * (hi[2] - lo[2]));
        cnt++;
    }
    return cnt ? sum / cnt : 0.0;
}

/* Tryptophan in water, the waters shuffled (whole molecules swap places) so
   their order says nothing of where they are: the update must regroup them,
   and the clustered licorice must give the depths and picks of a plain one. */
/* Point p's radius as md_render reads it: the float at radius + p * stride,
   or with a radius table the table's entry for the type in its low 16 bits */
static float* read_radius_table(fix_t* fx, const md_render_batch_t* b) {
    if (!b->points.radius_table) return NULL;
    return (float*)read_block(fx, b->points.radius_table, (uint64_t)b->points.radius_table_count * 4u + 4u);
}
static float radius_at(const float* r, uint32_t stride_floats, uint32_t p, const float* table, uint32_t count) {
    float v = r[(size_t)stride_floats * p];
    if (!table) return v;
    uint32_t w; memcpy(&w, &v, 4);
    return (w & 0xFFFFu) < count ? table[w & 0xFFFFu] : 0.0f;
}

/* The copy's radii follow their source across regroups: the types' (the
   regroup writes the small molecules' from their kinds) or the caller's
   array, switched with and without a regroup, no flags */
static int radii_wrong(fix_t* fx, const md_render_mol_topo_t* topo, md_gpu_addr_t dyn, const md_render_mol_atoms_t* atoms,
                       const md_render_mol_rep_t* rep, const float* expect) {
    md_render_batch_t b;
    if (md_render_mol_emit(topo, dyn, atoms, rep, &b, 1) != 1) return -1;
    uint8_t* blk = read_block(fx, dyn, topo->dynamic_size);
    uint8_t* ord = read_block(fx, topo->atom_order, (size_t)topo->point_count * 4);
    const float* r = (const float*)(blk + (b.points.radius - dyn));
    const uint32_t* order = (const uint32_t*)ord;
    float* table = read_radius_table(fx, &b);
    int wrong = 0;
    for (uint32_t p = 0; p < topo->point_count; ++p) wrong += radius_at(r, 4, p, table, b.points.radius_table_count) != expect[order[p]];
    free(blk); free(ord); free(table);
    return wrong;
}

UTEST(render, mol_radii_follow_their_source) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(load_system(&sys, &state, MD_UNITTEST_DATA_DIR "/tryptophan-md.gro"));
    const uint32_t n = (uint32_t)sys.atom.count;
    md_temp_scope_t temp = md_temp_begin();
    md_render_mol_plan_t plan;
    ASSERT_TRUE(md_render_mol_plan(&plan, &sys, NULL, temp));
    md_render_mol_topo_t topo;
    ASSERT_TRUE(md_render_mol_write(fx.gfx, &topo, &plan, md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, plan.size).gpu, NULL));
    md_temp_end(temp);
    md_render_mol_t mol = md_render_mol_create(fx.dev, NULL);
    ASSERT_TRUE(mol != NULL);
    md_gpu_addr_t dyn = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, topo.dynamic_size).gpu;
    md_render_mol_rep_t rep = {0};
    rep.key = 1; rep.type = MD_RENDER_MOL_REP_SPACEFILL; rep.params.spacefill.radius_scale = 1.0f;

    float* by_type = malloc(n * 4), *own = malloc(n * 4);
    for (uint32_t i = 0; i < n; ++i) {
        by_type[i] = sys.atom.type.radius[sys.atom.type_idx[i]];
        own[i] = 0.5f + 0.001f * (float)(i % 997);          /* no type's */
    }
    md_render_mol_atoms_t types = {0}, mine = {0};
    types.position = f_buf(&fx, state.xyz, n * 12);
    mine = types;
    mine.radius = f_buf(&fx, own, n * 4);
    const md_render_mol_update_flags_t RG = MD_RENDER_MOL_UPDATE_REGROUP;
    const struct { const md_render_mol_atoms_t* atoms; md_render_mol_update_flags_t flags; const float* expect; const char* what; } steps[] = {
        { &types, RG, by_type, "types, regroup" },
        { &types, 0,  by_type, "types" },
        { &mine,  0,  own,     "own, no regroup (the source changed)" },
        { &mine,  RG, own,     "own, regroup (the write's from the array)" },
        { &mine,  RG, own,     "own, regroup again" },
        { &types, 0,  by_type, "types, no regroup (the source changed back)" },
        { &types, RG, by_type, "types, regroup" },
    };
    for (size_t k = 0; k < sizeof(steps) / sizeof(steps[0]); ++k) {
        ASSERT_TRUE(md_render_mol_update(mol, fx.gfx, &topo, dyn, steps[k].atoms, steps[k].flags, NULL));
        const int wrong = radii_wrong(&fx, &topo, dyn, steps[k].atoms, &rep, steps[k].expect);
        if (wrong) printf("    %s: %d points with a wrong radius\n", steps[k].what, wrong);
        EXPECT_EQ(0, wrong);
    }
    /* A topology written again at the same address, of the same system
       version and the same number of points, but another one (the resident
       set differs by an atom): the copy made for the first is not the
       second's, nor are frames made for the first */
    {
        md_bitmask_t ra = md_bitmask_create(0, n, md_get_heap_allocator());
        md_bitmask_t rb = md_bitmask_create(0, n, md_get_heap_allocator());
        md_bitmask_set_range(&ra, 0, n); md_bitmask_clear_bit(&ra, 0);
        md_bitmask_set_range(&rb, 0, n); md_bitmask_clear_bit(&rb, 1);
        md_bitset_t sa = md_bitmask_view(&ra), sb = md_bitmask_view(&rb);
        md_temp_scope_t t2 = md_temp_begin();
        md_render_mol_plan_t pa, pb;
        ASSERT_TRUE(md_render_mol_plan(&pa, &sys, &sa, t2));
        ASSERT_TRUE(md_render_mol_plan(&pb, &sys, &sb, t2));
        ASSERT_EQ(pa.size, pb.size);
        ASSERT_EQ(pa.layout.point_count, pb.layout.point_count);
        md_gpu_addr_t base = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, pa.size).gpu;
        md_gpu_addr_t dyn2 = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, pa.layout.dynamic_size).gpu;
        md_render_mol_topo_t ta, tb;
        ASSERT_TRUE(md_render_mol_write(fx.gfx, &ta, &pa, base, NULL));
        ASSERT_TRUE(md_render_mol_update(mol, fx.gfx, &ta, dyn2, &types, 0, NULL));
        EXPECT_EQ(0, radii_wrong(&fx, &ta, dyn2, &types, &rep, by_type));
        md_render_mol_frames_t fr;
        md_gpu_addr_t frb = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, md_render_mol_frames_size(&ta, 1)).gpu;
        ASSERT_TRUE(md_render_mol_frames_init(&fr, &ta, frb, 1));
        {
            /* a frame in its slot, so only the topology can refuse it below */
            const float cell[9] = { 0 };
            uint16_t* q = malloc((size_t)n * 6);
            md_render_mol_frame_t f;
            md_render_mol_frame_quantize(&f, q, (const float*)state.xyz, 0, n, cell, &sa);
            ASSERT_TRUE(md_render_mol_frame_upload(mol, fx.gfx, &ta, dyn2, &fr, 0, &f, 0));
            md_render_mol_atoms_t af = {0};
            af.frames = &fr; af.frame_count = 1;
            EXPECT_TRUE(md_render_mol_update(mol, fx.gfx, &ta, dyn2, &af, 0, NULL));
            md_gpu_stream_sync(fx.gfx);
            free(q);
        }
        ASSERT_TRUE(md_render_mol_write(fx.gfx, &tb, &pb, base, NULL));
        EXPECT_NE(ta.instance, tb.instance);
        ASSERT_TRUE(md_render_mol_update(mol, fx.gfx, &tb, dyn2, &types, 0, NULL));
        const int wrong = radii_wrong(&fx, &tb, dyn2, &types, &rep, by_type);
        if (wrong) printf("    written again at the same base: %d points with a wrong radius\n", wrong);
        EXPECT_EQ(0, wrong);
        md_render_mol_atoms_t af = {0};
        af.frames = &fr; af.frame_count = 1;
        EXPECT_FALSE(md_render_mol_update(mol, fx.gfx, &tb, dyn2, &af, 0, NULL));   /* frames of the other one */
        md_temp_end(t2);
        md_gpu_free(fx.gfx, frb); md_gpu_free(fx.gfx, dyn2); md_gpu_free(fx.gfx, base);
        md_bitmask_free(&ra); md_bitmask_free(&rb);
    }
    free(by_type); free(own);
    md_gpu_free(fx.gfx, types.position); md_gpu_free(fx.gfx, mine.radius);
    md_gpu_free(fx.gfx, dyn); md_gpu_free(fx.gfx, topo.base);
    md_render_mol_destroy(mol);
    md_system_free(&sys);
    md_vm_arena_destroy(alloc);
    f_close(&fx);
}

UTEST(render, mol_clusters_follow_the_solvent) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(load_system(&sys, &state, MD_UNITTEST_DATA_DIR "/tryptophan-md.gro"));
    const uint32_t n = (uint32_t)sys.atom.count, nb = (uint32_t)(sys.bond.count + sys.bond.noncovalent_count);
    uint64_t content[MD_RENDER_CONTENT_COUNT] = {0};

    md_temp_scope_t temp = md_temp_begin();
    md_render_mol_plan_t plan;
    ASSERT_TRUE(md_render_mol_plan(&plan, &sys, NULL, temp));
    md_render_mol_topo_t topo;
    ASSERT_TRUE(md_render_mol_write(fx.gfx, &topo, &plan, md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, plan.size).gpu, content));
    md_temp_end(temp);
    EXPECT_EQ(27u, n - topo.small_atom_count);              /* the tryptophan */
    EXPECT_EQ(2156u, topo.small_molecule_count);            /* the waters */
    ASSERT_GT(topo.dynamic_size, 0u);

    md_render_mol_t mol = md_render_mol_create(fx.dev, NULL);
    ASSERT_TRUE(mol != NULL);
    md_gpu_addr_t dyn = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, topo.dynamic_size).gpu;

    /* Shuffle */
    float* xyz = malloc(n * 12);
    memcpy(xyz, state.xyz, n * 12);
    uint32_t seed = 7;
    const md_molecule_data_t* M = &sys.molecule;
    for (size_t k = M->count; k-- > 1;) {
        seed = seed * 1664525u + 1013904223u;
        size_t j = (seed >> 8) % (k + 1);
        uint32_t sk = M->offset[k + 1] - M->offset[k], sj = M->offset[j + 1] - M->offset[j];
        if (sk != sj || sk > MD_RENDER_MOL_SMALL) continue;
        for (uint32_t s = 0; s < sk; ++s) {
            float* a = xyz + 3 * M->atom_idx[M->offset[k] + s];
            float* b = xyz + 3 * M->atom_idx[M->offset[j] + s];
            for (int d = 0; d < 3; ++d) { float t = a[d]; a[d] = b[d]; b[d] = t; }
        }
    }
    md_gpu_addr_t P = f_buf(&fx, xyz, n * 12);
    content[MD_RENDER_CONTENT_GEOMETRY]++;
    md_render_mol_atoms_t atoms = {0};
    atoms.position = P;

    uint8_t* before = read_block(&fx, topo.base, topo.size);
    ASSERT_TRUE(md_render_mol_update(mol, fx.gfx, &topo, dyn, &atoms, MD_RENDER_MOL_UPDATE_REGROUP, content));
    uint8_t* after = read_block(&fx, topo.base, topo.size);
    /* the block's own image: rebase the addresses to offsets */
    md_render_mol_topo_t rel = topo;
    rel.atom_order -= topo.base; rel.mixed_point -= topo.base; rel.atom_cluster_offset -= topo.base; rel.bond_pairs -= topo.base; rel.bond_ids -= topo.base;
    for (int c = 0; c < MD_RENDER_MOL_BOND_CLASS_COUNT; ++c) rel.bond_cluster_offset[c] -= topo.base;
    EXPECT_EQ(0, check_topo_image(&sys, &rel, after, NULL));
    double e0 = grid_cluster_extent(&rel, before, xyz), e1 = grid_cluster_extent(&rel, after, xyz);
    printf("    grid clusters: mean extent %.1f A in molecule order, %.1f A after the update\n", e0, e1);
    EXPECT_LT(e1, 0.4 * e0);
    /* again from the regrouped order: the same cells (the order within one is the atomics') */
    ASSERT_TRUE(md_render_mol_update(mol, fx.gfx, &topo, dyn, &atoms, MD_RENDER_MOL_UPDATE_REGROUP, content));
    uint8_t* again = read_block(&fx, topo.base, topo.size);
    EXPECT_EQ(0, check_topo_image(&sys, &rel, again, NULL));
    EXPECT_EQ(0, memcmp(after + rel.atom_cluster_offset, again + rel.atom_cluster_offset, (topo.atom_cluster_count + 1) * 4));
    /* Velocities: copied into the point order when given (a per-frame block
       of dynamic_size_velocity), and emit hands md_render that copy */
    {
        EXPECT_GT(topo.dynamic_size_velocity, topo.dynamic_size);
        md_gpu_addr_t dyn_v = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, topo.dynamic_size_velocity).gpu;
        uint32_t* vel = malloc((size_t)n * 8);
        for (uint32_t i = 0; i < n; ++i) { vel[2 * i] = 0x10000u + i; vel[2 * i + 1] = 0x20000u + i; }
        md_render_mol_atoms_t av = atoms;
        av.velocity = f_buf(&fx, vel, (size_t)n * 8);
        /* no regroup: that would reorder the topology under `dyn`, whose copy is drawn below */
        ASSERT_TRUE(md_render_mol_update(mol, fx.gfx, &topo, dyn_v, &av, 0, content));
        md_render_mol_rep_t rv = {0};
        rv.key = 9; rv.type = MD_RENDER_MOL_REP_SPACEFILL; rv.params.spacefill.radius_scale = 1.0f;
        md_render_batch_t bv;
        ASSERT_EQ(1u, md_render_mol_emit(&topo, dyn_v, &av, &rv, &bv, 1));
        ASSERT_TRUE(bv.points.velocity != 0);
        uint8_t* blk  = read_block(&fx, dyn_v, topo.dynamic_size_velocity);
        uint8_t* topb = read_block(&fx, topo.base, topo.size);
        const uint32_t* cv    = (const uint32_t*)(blk + (bv.points.velocity - dyn_v));
        const uint32_t* order = (const uint32_t*)(topb + (topo.atom_order - topo.base));
        int wrong = 0;
        for (uint32_t p = 0; p < n; ++p) wrong += cv[2 * p] != vel[2 * order[p]] || cv[2 * p + 1] != vel[2 * order[p] + 1];
        EXPECT_EQ(0, wrong);
        /* and without velocities, no velocity for md_render */
        ASSERT_EQ(1u, md_render_mol_emit(&topo, dyn_v, &atoms, &rv, &bv, 1));
        EXPECT_TRUE(bv.points.velocity == 0);
        free(blk); free(topb); free(vel);
        md_gpu_free(fx.gfx, av.velocity);
        md_gpu_free(fx.gfx, dyn_v);
    }

    /* the radii of the particle types */
    {
        const float* r = (const float*)(again + (topo.type_radius - topo.base));
        int wrong = 0;
        for (size_t t = 0; t < sys.atom.type.count; ++t) wrong += r[t] != sys.atom.type.radius[t];
        EXPECT_EQ(0, wrong);
        EXPECT_EQ((uint32_t)sys.atom.type.count, topo.type_count);
    }

    /* Licorice: clustered (from emit) against plain batches in system order */
    md_render_mol_rep_t rep = {0};
    rep.key  = 7;
    rep.type = MD_RENDER_MOL_REP_LICORICE;
    rep.params.licorice.radius = 0.3f;
    rep.bond_classes = MD_RENDER_MOL_BONDS_ALL;
    md_render_batch_t batches[8];
    uint32_t nbat = md_render_mol_emit(&topo, dyn, &atoms, &rep, batches, 8);
    ASSERT_GE(nbat, 2u);
    EXPECT_EQ(topo.atom_cluster_count, batches[nbat - 1].cluster_count);

    uint32_t* plain_pairs = malloc(nb * 8);
    for (uint32_t i = 0; i < nb; ++i) { plain_pairs[2 * i] = (uint32_t)sys.bond.pairs[i].idx[0]; plain_pairs[2 * i + 1] = (uint32_t)sys.bond.pairs[i].idx[1]; }
    md_render_batch_t plain[2];
    memset(plain, 0, sizeof(plain));
    plain[0].key = 1; plain[0].kind = MD_RENDER_SPHERES; plain[0].count = n;
    plain[0].points.count = n; plain[0].points.position = P;
    plain[0].style.spheres.radius_add = 0.3f;
    plain[1].key = 2; plain[1].kind = MD_RENDER_CAPSULES; plain[1].count = nb;
    plain[1].points.count = n; plain[1].points.position = P; plain[1].index = f_buf(&fx, plain_pairs, nb * 8);
    plain[1].style.capsules.radius = 0.3f;

    /* the 40 A box from outside, in perspective */
    const uint32_t W = 160, H = 120;
    cam_t c = { .tan_half = 0.5f, .aspect = (float)W / H, .n = 1.0f };
    md_render_camera_t cam = make_camera(&c);
    cam.view[12] = -20.0f; cam.view[13] = -20.0f; cam.view[14] = -85.0f;
    md_render_view_t ref = md_render_view_create(fx.r, "reference");
    md_render_frame_t fr = frame_of(W, H, cam, plain, 2);
    fr.options = MD_RENDER_NO_OCCLUSION_CULLING | MD_RENDER_NO_COMPUTE_RASTER;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, ref, &fr, NULL));
    md_render_frame_t f = frame_of(W, H, cam, batches, nbat);
    memcpy(f.content, content, sizeof(content));
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    f.content[MD_RENDER_CONTENT_GEOMETRY]++;                 /* again, culling against the pyramid */
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    md_gpu_stream_sync(fx.gfx);
    md_render_stats_t st;
    ASSERT_TRUE(md_render_view_stats(fx.view, &st));
    printf("    clusters %u, expanded %u + %u, occluded %u; elements %llu, occluded %u\n", st.clusters,
           st.clusters_expanded[0], st.clusters_expanded[1], st.clusters_occluded, (unsigned long long)st.elements, st.occluded);
    int covered = 0;
    EXPECT_EQ(0, compare_views_ex(&fx, fx.view, ref, W, H, &covered, true));
    EXPECT_GT(covered, (int)(W * H / 4));

    /* Picks: the same atom or bond under the same pixel */
    int picked = 0, agree = 0;
    for (uint32_t y = 5; y < H; y += 23) for (uint32_t x = 7; x < W; x += 29) {
        md_render_pick_t pa, pb;
        ASSERT_TRUE(f_pick(&fx, x, y, &pa));
        if (!md_render_view_pick(fx.gfx, ref, x, y)) continue;
        md_gpu_stream_sync(fx.gfx);
        ASSERT_TRUE(md_render_view_pick_result(ref, &pb));
        if (pb.key == 0) { EXPECT_EQ(0u, (unsigned)pa.key); continue; }
        md_render_mol_hit_t hit;
        ASSERT_TRUE(md_render_mol_decode_pick(&pa, &hit));
        picked++;
        if (pb.kind == MD_RENDER_SPHERES) agree += hit.kind == MD_RENDER_MOL_HIT_ATOM && hit.atom == pb.attr;
        else                              agree += hit.kind == MD_RENDER_MOL_HIT_BOND && hit.bond == pb.id;
    }
    EXPECT_GT(picked, 10);
    EXPECT_EQ(picked, agree);

    /* Spacefill: the radii of the particle types, then of an edited table */
    float* rad = malloc(n * 4);
    for (int pass = 0; pass < 2; ++pass) {
        const float scale = pass ? 0.5f : 1.0f;
        if (pass) {
            float* table = malloc(sys.atom.type.count * 4);
            for (size_t t = 0; t < sys.atom.type.count; ++t) table[t] = sys.atom.type.radius[t] * scale;
            atoms.radius_table = f_buf(&fx, table, sys.atom.type.count * 4);
            free(table);
            ASSERT_TRUE(md_render_mol_update(mol, fx.gfx, &topo, dyn, &atoms, MD_RENDER_MOL_UPDATE_RADII, content));   /* refills the copy too */
        }
        /* a plain update (positions only) keeps the copy's radii: its w, not the source */
        ASSERT_TRUE(md_render_mol_update(mol, fx.gfx, &topo, dyn, &atoms, 0, content));
        for (uint32_t i = 0; i < n; ++i) rad[i] = sys.atom.type.radius[sys.atom.type_idx[i]] * scale;
        md_gpu_addr_t R = f_buf(&fx, rad, n * 4);
        rep.type = MD_RENDER_MOL_REP_SPACEFILL;
        rep.params.spacefill.radius_scale = 1.0f;
        nbat = md_render_mol_emit(&topo, dyn, &atoms, &rep, batches, 8);
        ASSERT_EQ(1u, nbat);
        md_render_batch_t ps = plain[0];
        ps.points.radius = R;
        ps.style.spheres.radius_add = 0.0f;
        ps.style.spheres.radius_scale = 1.0f;
        fr.batches = &ps; fr.batch_count = 1;
        fr.content[MD_RENDER_CONTENT_GEOMETRY]++;
        ASSERT_TRUE(md_render_view_draw(fx.gfx, ref, &fr, NULL));
        f.batches = batches; f.batch_count = nbat;
        memcpy(f.content, content, sizeof(content));
        f.content[MD_RENDER_CONTENT_GEOMETRY] += 1000;
        ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
        EXPECT_EQ(0, compare_views_ex(&fx, fx.view, ref, W, H, &covered, true));
        EXPECT_GT(covered, (int)(W * H / 4));
        md_gpu_stream_sync(fx.gfx);
        md_gpu_free(fx.gfx, R);
    }
    free(rad);
    md_gpu_free(fx.gfx, atoms.radius_table);

    md_render_view_destroy(fx.gfx, ref);
    free(xyz); free(before); free(after); free(again); free(plain_pairs);
    md_gpu_free(fx.gfx, plain[1].index);
    md_gpu_free(fx.gfx, P);
    md_gpu_free(fx.gfx, dyn);
    md_gpu_free(fx.gfx, topo.base);
    md_render_mol_destroy(mol);
    md_system_free(&sys);
    f_close(&fx);
}

/* ---- Frames: resident on the GPU, interpolated there ---------------------------- */

/* All atoms, or (partial) the bonded structures, every third small molecule
   and the first atom of every fifth: whole molecules, parts of molecules and
   molecules left out */
static md_bitmask_t make_resident(const md_system_t* sys, bool partial) {
    const uint32_t n = (uint32_t)sys->atom.count;
    md_bitmask_t mask = md_bitmask_create(0, n, md_get_heap_allocator());
    const md_molecule_data_t* M = &sys->molecule;
    uint32_t w = 0;
    for (size_t m = 0; m < M->count; ++m) {
        uint32_t b = M->offset[m], e = M->offset[m + 1];
        bool small = e - b <= MD_RENDER_MOL_SMALL;
        if (small) w++;
        for (uint32_t k = b; k < e; ++k) {
            bool keep = !partial || !small || w % 3 == 0 || (w % 5 == 0 && k == b);
            if (keep) md_bitmask_set_bit(&mask, (uint32_t)M->atom_idx[k]);
        }
    }
    return mask;
}


static float half_to_float(uint32_t h) {
    uint32_t sgn = (h >> 15) & 1u, e = (h >> 10) & 31u, m = h & 1023u;
    float v = e == 0 ? ldexpf((float)m, -24) : e == 31 ? (m ? NAN : INFINITY) : ldexpf((float)(m | 1024u), (int)e - 25);
    return sgn ? -v : v;
}

/* The tryptophan box over five frames: every atom on its own curve, wrapped
   into the box, so many cross it between frames. Uploads that regroup (one
   into a slot that holds a frame) must leave every resident frame in the new
   point order, and nearest, linear and cubic must give the CPU's answer from
   the same quantized frames, minimum image included. `partial`: only some
   atoms resident (the residue, every third water, the oxygen of every
   fifth): only those are quantized, uploaded, regrouped and drawn. */
/* in md_render_mol.c, not in its header: the frames' gather from at most `groups` groups */
void md_render_mol_debug_frame_groups(uint32_t groups);

static void frames_case(int* utest_result, bool partial) {
    fix_t fx;
    if (!f_open(&fx)) { printf("    skipped: %s\n", f_skip_reason()); return; }
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(load_system(&sys, &state, MD_UNITTEST_DATA_DIR "/tryptophan-md.gro"));
    const uint32_t n = (uint32_t)sys.atom.count;
    ASSERT_TRUE(md_unitcell_is_orthorhombic(&state.unitcell));
    const float L[3] = { (float)state.unitcell.x, (float)state.unitcell.y, (float)state.unitcell.z };
    float cell[9];
    md_render_mol_frame_cell(cell, &state.unitcell);
    EXPECT_EQ(L[0], cell[0]); EXPECT_EQ(L[1], cell[4]); EXPECT_EQ(L[2], cell[8]);

    /* the resident atoms */
    md_bitmask_t mask = make_resident(&sys, partial);
    md_bitset_t rset = md_bitmask_view(&mask);
    const md_bitset_t* resident = partial ? &rset : NULL;
    uint32_t* rlist = malloc((size_t)n * 4);           /* the resident atoms in order: a frame's rows */
    uint32_t nres = (uint32_t)md_bitset_extract_indices(rlist, n, rset);
    printf("    %s: %u of %u atoms resident\n", partial ? "partial" : "all", nres, n);

    md_temp_scope_t temp = md_temp_begin();
    md_render_mol_plan_t plan;
    ASSERT_TRUE(md_render_mol_plan(&plan, &sys, resident, temp));
    EXPECT_EQ(nres, plan.layout.point_count);
    EXPECT_EQ(0, check_topo_image(&sys, &plan.layout, plan.data, resident));
    md_render_mol_topo_t topo;
    ASSERT_TRUE(md_render_mol_write(fx.gfx, &topo, &plan, md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, plan.size).gpu, NULL));
    md_temp_end(temp);
    md_render_mol_t mol = md_render_mol_create(fx.dev, NULL);
    ASSERT_TRUE(mol != NULL);
    md_gpu_addr_t dyn = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, topo.dynamic_size_velocity).gpu;

    /* Five frames: x(k) = x + v k + c k^2, wrapped; the waters shuffled
       first, so their file order says nothing of where they are */
    enum { NF = 5 };
    float* xyz = malloc((size_t)n * 12);
    memcpy(xyz, state.xyz, (size_t)n * 12);
    uint32_t seed = 11;
    #define RND() (((seed = seed * 1664525u + 1013904223u) >> 8) / 16777216.0f)
    const md_molecule_data_t* M = &sys.molecule;
    for (size_t k = M->count; k-- > 1;) {
        size_t j = (size_t)(RND() * (k + 1));
        uint32_t sk = M->offset[k + 1] - M->offset[k], sj = M->offset[j + 1] - M->offset[j];
        if (sk != sj || sk > MD_RENDER_MOL_SMALL) continue;
        for (uint32_t q = 0; q < sk; ++q) {
            float* a = xyz + 3 * M->atom_idx[M->offset[k] + q];
            float* b = xyz + 3 * M->atom_idx[M->offset[j] + q];
            for (int d = 0; d < 3; ++d) { float t = a[d]; a[d] = b[d]; b[d] = t; }
        }
    }
    float* fpos[NF];
    uint16_t* qs[NF];
    md_render_mol_frame_t frames[NF];
    float* deq[NF];                                /* what the GPU holds, atom order */
    float* v = malloc((size_t)n * 24);
    for (uint32_t i = 0; i < 6 * n; ++i) v[i] = (RND() - 0.5f) * (i % 2 ? 0.3f : 1.6f);
    int crossings = 0, quant_bad = 0;
    for (int f = 0; f < NF; ++f) {
        fpos[f] = malloc((size_t)n * 12);
        deq[f]  = malloc((size_t)n * 12);
        for (uint32_t i = 0; i < n; ++i) for (int d = 0; d < 3; ++d) {
            float x = xyz[3 * i + d] + v[6 * i + 2 * d] * f + v[6 * i + 2 * d + 1] * f * f;
            fpos[f][3 * i + d] = x - L[d] * floorf(x / L[d]);
            if (f && fabsf(fpos[f][3 * i + d] - fpos[f - 1][3 * i + d]) > 0.5f * L[d]) crossings++;
        }
        qs[f] = malloc((size_t)n * 6);
        md_render_mol_frame_quantize(&frames[f], qs[f], fpos[f], 0, n, cell, resident);
        for (uint32_t r = 0; r < nres; ++r) for (int d = 0; d < 3; ++d) {
            const uint32_t i = rlist[r];
            deq[f][3 * i + d] = frames[f].origin[d] + (float)qs[f][3 * r + d] * frames[f].scale[d];
            quant_bad += fabsf(deq[f][3 * i + d] - fpos[f][3 * i + d]) > 0.5f * frames[f].scale[d] + 1e-5f;
        }
    }
    printf("    %d crossings of the box between frames; step %.4f A\n", crossings, frames[0].scale[0]);
    EXPECT_GT(crossings, 100);
    EXPECT_EQ(0, quant_bad);
    /* the vector path and the scalar one agree: the same frame through a
       stride the vector path does not take */
    if (!partial) {
        float* wide = malloc((size_t)n * 16);
        uint16_t* q2 = malloc((size_t)n * 6);
        for (uint32_t i = 0; i < n; ++i) { memcpy(wide + 4 * i, fpos[2] + 3 * i, 12); wide[4 * i + 3] = NAN; }
        md_render_mol_frame_t f2;
        md_render_mol_frame_quantize(&f2, q2, wide, 16, n, cell, NULL);
        EXPECT_EQ(0, memcmp(f2.origin, frames[2].origin, sizeof(f2.origin)));
        EXPECT_EQ(0, memcmp(f2.scale, frames[2].scale, sizeof(f2.scale)));
        int diff = 0;
        for (uint32_t k = 0; k < 3 * n; ++k) diff += abs((int)q2[k] - (int)qs[2][k]) > 1;
        EXPECT_EQ(0, diff);
        free(wide); free(q2);
    }
    /* frames 3 and 4 from GPU-visible memory the host quantized into (the
       unified-memory path: nothing copied), the others from the host */
    md_gpu_mem_t up[2];
    for (int k = 0; k < 2; ++k) {
        up[k] = md_gpu_malloc(fx.gfx, MD_GPU_MEM_HOST_WRITE, md_render_mol_frame_size(&topo));
        ASSERT_TRUE(up[k].cpu != NULL);
        md_render_mol_frame_quantize(&frames[3 + k], up[k].cpu, fpos[3 + k], 0, n, cell, resident);
        frames[3 + k].q_gpu = up[k].gpu;
    }

    md_render_mol_frames_t fr;
    const uint64_t fsize = md_render_mol_frames_size(&topo, 4);
    ASSERT_GE(fsize, (uint64_t)nres * 24);           /* four slots of 6 bytes per atom, no staging */
    ASSERT_LT(fsize, (uint64_t)nres * 29 + 4096);
    ASSERT_TRUE(md_render_mol_frames_init(&fr, &topo, md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, fsize).gpu, 4));
    const md_render_mol_update_flags_t RG = MD_RENDER_MOL_UPDATE_REGROUP;
    /* frames 0..3 into slots 0..3, regrouping with the first and the last */
    uint8_t* order0 = read_block(&fx, topo.atom_order, (size_t)nres * 4);
    ASSERT_TRUE(md_render_mol_frame_upload(mol, fx.gfx, &topo, dyn, &fr, 0, &frames[0], RG));
    uint8_t* order1 = read_block(&fx, topo.atom_order, (size_t)nres * 4);
    ASSERT_TRUE(md_render_mol_frame_upload(mol, fx.gfx, &topo, 0,   &fr, 1, &frames[1], 0));
    ASSERT_TRUE(md_render_mol_frame_upload(mol, fx.gfx, &topo, 0,   &fr, 2, &frames[2], 0));
    ASSERT_TRUE(md_render_mol_frame_upload(mol, fx.gfx, &topo, dyn, &fr, 3, &frames[3], RG));
    uint8_t* order2 = read_block(&fx, topo.atom_order, (size_t)nres * 4);
    EXPECT_NE(0, memcmp(order0, order1, (size_t)nres * 4));
    EXPECT_NE(0, memcmp(order1, order2, (size_t)nres * 4));
    EXPECT_EQ(15u, fr.valid);

    /* the slots' frames, by slot: which of ours each holds */
    int held[4] = { 0, 1, 2, 3 };
    md_render_mol_atoms_t atoms = {0};
    atoms.frames = &fr;
    md_render_mol_rep_t rep = {0};
    rep.key = 3; rep.type = MD_RENDER_MOL_REP_SPACEFILL; rep.params.spacefill.radius_scale = 1.0f;

    /* nearest, linear, cubic; then frame 4 into slot 0 (which held frame 0),
       regrouping, and cubic over the four latest */
    const struct { uint32_t count; uint32_t slot[4]; float t, dt; bool replace; } cases[] = {
        { 1, { 2 },          0.0f, 0.0f, false },
        { 2, { 1, 2 },       0.6f, 2.0f, false },
        { 4, { 0, 1, 2, 3 }, 0.3f, 2.0f, false },
        { 4, { 1, 2, 3, 0 }, 0.8f, 0.5f, true  },
    };
    for (size_t ci = 0; ci < sizeof(cases) / sizeof(cases[0]); ++ci) {
        if (cases[ci].replace) {
            ASSERT_TRUE(md_render_mol_frame_upload(mol, fx.gfx, &topo, dyn, &fr, 0, &frames[4], RG));
            held[0] = 4;
        }
        atoms.frame_count = cases[ci].count;
        memcpy(atoms.frame_slot, cases[ci].slot, sizeof(atoms.frame_slot));
        md_render_mol_frame_weights(atoms.frame_w, atoms.frame_dw, cases[ci].count, cases[ci].t, cases[ci].dt, 0.5f);
        ASSERT_TRUE(md_render_mol_update(mol, fx.gfx, &topo, dyn, &atoms, 0, NULL));
        md_render_batch_t b;
        ASSERT_EQ(1u, md_render_mol_emit(&topo, dyn, &atoms, &rep, &b, 1));
        const bool with_v = cases[ci].count > 1;
        EXPECT_EQ(with_v, b.points.velocity != 0);

        uint8_t* blk = read_block(&fx, dyn, topo.dynamic_size_velocity);
        uint8_t* ord = read_block(&fx, topo.atom_order, (size_t)nres * 4);
        const float* pts = (const float*)(blk + (b.points.position - dyn));
        const float* rad = (const float*)(blk + (b.points.radius - dyn));
        float* rtab = read_radius_table(&fx, &b);
        EXPECT_EQ(16u, b.points.stride);
        EXPECT_EQ(16u, b.points.radius_stride);
        EXPECT_TRUE(rtab != NULL);                  /* the types' table: w is the type */
        const uint32_t* vel = with_v ? (const uint32_t*)(blk + (b.points.velocity - dyn)) : NULL;
        const uint32_t* order = (const uint32_t*)ord;
        double pos_err = 0, vel_err = 0, vmax = 0;
        int rad_bad = 0;
        for (uint32_t p = 0; p < nres; ++p) {
            const uint32_t i = order[p];
            double x[4][3];
            for (uint32_t k = 0; k < cases[ci].count; ++k)
                for (int d = 0; d < 3; ++d) x[k][d] = deq[held[cases[ci].slot[k]]][3 * i + d];
            /* the minimum image, as documented: of [0] for two, of [1] (and [3] of [2]) for four */
            for (int d = 0; d < 3; ++d) {
                #define MI(a, ref) (a) = (ref) + ((a) - (ref)) - L[d] * roundf((float)(((a) - (ref)) / L[d]))
                if (cases[ci].count == 2) { MI(x[1][d], x[0][d]); }
                if (cases[ci].count == 4) { MI(x[0][d], x[1][d]); MI(x[2][d], x[1][d]); MI(x[3][d], x[2][d]); }
                #undef MI
                double pr = 0, vr = 0;
                for (uint32_t k = 0; k < cases[ci].count; ++k) { pr += atoms.frame_w[k] * x[k][d]; vr += atoms.frame_dw[k] * x[k][d]; }
                double e = fabs(pts[4 * p + d] - pr);
                if (e > pos_err) pos_err = e;
                if (with_v) {
                    uint32_t h = d == 0 ? vel[2 * p] & 0xFFFFu : d == 1 ? vel[2 * p] >> 16 : vel[2 * p + 1] & 0xFFFFu;
                    /* half precision, and the float sum of terms of 45 A: relative and absolute */
                    double ev = fabs(half_to_float(h) - vr) / (2e-3 * fabs(vr) + 1e-4);
                    if (ev > vel_err) vel_err = ev;
                    if (fabs(vr) > vmax) vmax = fabs(vr);
                }
            }
            rad_bad += radius_at(rad, 4, p, rtab, b.points.radius_table_count) != sys.atom.type.radius[sys.atom.type_idx[i]];
        }
        printf("    %u frame(s), t %.1f: position error %.2g A, velocity error %.2g of the tolerance (max |v| %.2f)\n",
               cases[ci].count, cases[ci].t, pos_err, vel_err, vmax);
        EXPECT_LT(pos_err, 1e-3);
        if (with_v) EXPECT_LT(vel_err, 1.0);
        EXPECT_EQ(0, rad_bad);
        {
            /* the same from one group, its threads striding over every point */
            md_render_mol_debug_frame_groups(1);
            ASSERT_TRUE(md_render_mol_update(mol, fx.gfx, &topo, dyn, &atoms, 0, NULL));
            md_render_mol_debug_frame_groups(0);
            uint8_t* blk1 = read_block(&fx, dyn, topo.dynamic_size_velocity);
            const size_t po = (size_t)(b.points.position - dyn);
            EXPECT_EQ(0, memcmp(blk + po, blk1 + po, (size_t)nres * 16));
            if (with_v) { const size_t vo = (size_t)(b.points.velocity - dyn); EXPECT_EQ(0, memcmp(blk + vo, blk1 + vo, (size_t)nres * 8)); }
            free(blk1);
        }
        free(blk); free(ord); free(rtab);
    }
    /* the topology after the regroups: whole and consistent */
    uint8_t* img = read_block(&fx, topo.base, topo.size);
    md_render_mol_topo_t rel = topo;
    rel.atom_order -= topo.base; rel.mixed_point -= topo.base; rel.atom_cluster_offset -= topo.base; rel.bond_pairs -= topo.base; rel.bond_ids -= topo.base;
    for (int c = 0; c < MD_RENDER_MOL_BOND_CLASS_COUNT; ++c) rel.bond_cluster_offset[c] -= topo.base;
    EXPECT_EQ(0, check_topo_image(&sys, &rel, img, resident));
    /* refused: a regroup in the update, a slot without a frame, a cell with a
       not along x (the minimum image is GROMACS's, for that form) */
    EXPECT_FALSE(md_render_mol_update(mol, fx.gfx, &topo, dyn, &atoms, RG, NULL));
    {
        md_render_mol_frame_t tilted = frames[0];
        tilted.cell[1] = 0.5f;
        EXPECT_FALSE(md_render_mol_frame_upload(mol, fx.gfx, &topo, 0, &fr, 0, &tilted, 0));
    }
    md_render_mol_frames_t empty;
    ASSERT_TRUE(md_render_mol_frames_init(&empty, &topo, fr.base, 4));
    atoms.frames = &empty;
    EXPECT_FALSE(md_render_mol_update(mol, fx.gfx, &topo, dyn, &atoms, 0, NULL));
    #undef RND

    md_gpu_stream_sync(fx.gfx);
    free(img); free(order0); free(order1); free(order2); free(v); free(xyz);
    for (int f = 0; f < NF; ++f) { free(fpos[f]); free(deq[f]); free(qs[f]); }
    md_gpu_free(fx.gfx, up[0].gpu);
    md_gpu_free(fx.gfx, up[1].gpu);
    md_gpu_free(fx.gfx, fr.base);
    md_gpu_free(fx.gfx, dyn);
    md_gpu_free(fx.gfx, topo.base);
    md_render_mol_destroy(mol);
    md_system_free(&sys);
    md_bitmask_free(&mask);
    free(rlist);
    f_close(&fx);
}

UTEST(render, mol_frames_interpolate) {
    frames_case(utest_result, false);
}

UTEST(render, mol_frames_resident_set) {
    frames_case(utest_result, true);
}

/* Leaving atoms out of the topology draws what hiding them draws: licorice
   from a topology of the resident atoms against plain batches of every atom
   with the others hidden (and the bonds between resident atoms), depth by
   depth, and picks of the same atoms and bonds */
UTEST(render, mol_resident_set_draws_what_hiding_draws) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(load_system(&sys, &state, MD_UNITTEST_DATA_DIR "/tryptophan-md.gro"));
    const uint32_t n = (uint32_t)sys.atom.count, nb = (uint32_t)sys.bond.count;
    md_bitmask_t mask = make_resident(&sys, true);
    md_bitset_t rset = md_bitmask_view(&mask);
    const uint32_t nres = (uint32_t)md_bitset_count(rset);

    md_temp_scope_t temp = md_temp_begin();
    md_render_mol_plan_t plan;
    ASSERT_TRUE(md_render_mol_plan(&plan, &sys, &rset, temp));
    md_render_mol_topo_t topo;
    ASSERT_TRUE(md_render_mol_write(fx.gfx, &topo, &plan, md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, plan.size).gpu, NULL));
    md_temp_end(temp);
    EXPECT_EQ(nres, topo.point_count);
    md_render_mol_t mol = md_render_mol_create(fx.dev, NULL);
    ASSERT_TRUE(mol != NULL);
    md_gpu_addr_t dyn = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, topo.dynamic_size).gpu;
    md_gpu_addr_t P = f_buf(&fx, state.xyz, (size_t)n * 12);
    md_render_mol_atoms_t atoms = {0};
    atoms.position = P;
    ASSERT_TRUE(md_render_mol_update(mol, fx.gfx, &topo, dyn, &atoms, MD_RENDER_MOL_UPDATE_REGROUP, NULL));

    md_render_mol_rep_t rep = {0};
    rep.key  = 7;
    rep.type = MD_RENDER_MOL_REP_LICORICE;
    rep.params.licorice.radius = 0.3f;
    rep.bond_classes = MD_RENDER_MOL_BONDS_ALL;
    md_render_batch_t batches[8];
    uint32_t nbat = md_render_mol_emit(&topo, dyn, &atoms, &rep, batches, 8);
    ASSERT_GE(nbat, 2u);
    EXPECT_EQ(nres, batches[nbat - 1].count);
    EXPECT_EQ(n, batches[nbat - 1].points.attr_count);

    /* the reference: every atom, the others hidden; the bonds between resident atoms */
    uint32_t* bits = calloc((n + 31) / 32, 4);
    for (uint32_t i = 0; i < n; ++i) if (md_bitset_test(rset, i)) bits[i >> 5] |= 1u << (i & 31);
    uint32_t* pairs = malloc((size_t)nb * 8);
    uint32_t* ids   = malloc((size_t)nb * 4);
    uint32_t kept = 0;
    for (uint32_t i = 0; i < nb; ++i) {
        uint32_t a = (uint32_t)sys.bond.pairs[i].idx[0], b = (uint32_t)sys.bond.pairs[i].idx[1];
        if (!md_bitset_test(rset, a) || !md_bitset_test(rset, b)) continue;
        pairs[2 * kept] = a; pairs[2 * kept + 1] = b; ids[kept] = i; kept++;
    }
    EXPECT_EQ(kept, topo.bond_count);
    md_render_batch_t plain[2];
    memset(plain, 0, sizeof(plain));
    plain[0].key = 1; plain[0].kind = MD_RENDER_SPHERES; plain[0].count = n;
    plain[0].points.count = n; plain[0].points.position = P;
    plain[0].points.visible[0] = f_buf(&fx, bits, (n + 31) / 32 * 4);
    plain[0].style.spheres.radius_add = 0.3f;
    plain[1].key = 2; plain[1].kind = MD_RENDER_CAPSULES; plain[1].count = kept;
    plain[1].points.count = n; plain[1].points.position = P;
    plain[1].index = f_buf(&fx, pairs, (size_t)kept * 8);
    plain[1].id    = f_buf(&fx, ids, (size_t)kept * 4);
    plain[1].style.capsules.radius = 0.3f;

    const uint32_t W = 160, H = 120;
    cam_t c = { .tan_half = 0.5f, .aspect = (float)W / H, .n = 1.0f };
    md_render_camera_t cam = make_camera(&c);
    cam.view[12] = -20.0f; cam.view[13] = -20.0f; cam.view[14] = -85.0f;
    md_render_view_t ref = md_render_view_create(fx.r, "reference");
    md_render_frame_t fr = frame_of(W, H, cam, plain, 2);
    fr.options = MD_RENDER_NO_OCCLUSION_CULLING | MD_RENDER_NO_COMPUTE_RASTER;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, ref, &fr, NULL));
    md_render_frame_t f = frame_of(W, H, cam, batches, nbat);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    f.content[MD_RENDER_CONTENT_GEOMETRY]++;                 /* again, culling against the pyramid */
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    md_gpu_stream_sync(fx.gfx);
    int covered = 0;
    EXPECT_EQ(0, compare_views_ex(&fx, fx.view, ref, W, H, &covered, true));
    printf("    %u of %u atoms resident, %u of %u bonds; %d pixels covered\n", nres, n, kept, nb, covered);
    EXPECT_GT(covered, (int)(W * H / 8));

    int picked = 0, agree = 0;
    for (uint32_t y = 3; y < H; y += 11) for (uint32_t x = 5; x < W; x += 13) {
        md_render_pick_t pa, pb;
        ASSERT_TRUE(f_pick(&fx, x, y, &pa));
        if (!md_render_view_pick(fx.gfx, ref, x, y)) continue;
        md_gpu_stream_sync(fx.gfx);
        ASSERT_TRUE(md_render_view_pick_result(ref, &pb));
        if (pb.key == 0) { EXPECT_EQ(0u, (unsigned)pa.key); continue; }
        md_render_mol_hit_t hit;
        ASSERT_TRUE(md_render_mol_decode_pick(&pa, &hit));
        picked++;
        if (pb.kind == MD_RENDER_SPHERES) agree += hit.kind == MD_RENDER_MOL_HIT_ATOM && hit.atom == pb.attr;
        else                              agree += hit.kind == MD_RENDER_MOL_HIT_BOND && hit.bond == pb.id;
    }
    EXPECT_GT(picked, 10);
    EXPECT_EQ(picked, agree);

    md_gpu_stream_sync(fx.gfx);
    md_render_view_destroy(fx.gfx, ref);
    free(bits); free(pairs); free(ids);
    md_gpu_free(fx.gfx, dyn);
    md_gpu_free(fx.gfx, topo.base);
    md_render_mol_destroy(mol);
    md_system_free(&sys);
    md_bitmask_free(&mask);
    f_close(&fx);
}

/* ---- The cartoon ------------------------------------------------------------------
   The protein backbone through md_render's splines: control points at the CAs
   (the topology), axes oriented by the update with memory across updates,
   profiles from the secondary structure; drawn and picked as residues. */

static uint32_t lcg_next(uint32_t* s) { *s = *s * 1664525u + 1013904223u; return *s; }
static float lcg_gauss(uint32_t* s) {
    float u1 = ((lcg_next(s) >> 8) + 1u) * (1.0f / 16777217.0f), u2 = (lcg_next(s) >> 8) * (1.0f / 16777216.0f);
    return sqrtf(-2.0f * logf(u1)) * cosf(6.2831853f * u2);
}

UTEST(render, mol_cartoon_backbone) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(load_system(&sys, &state, MD_UNITTEST_DATA_DIR "/tubulin-A-B.pdb"));
    const md_protein_backbone_data_t* PB = &sys.protein_backbone;
    ASSERT_GT((uint32_t)PB->segment.count, 50u);
    uint64_t content[MD_RENDER_CONTENT_COUNT] = {0};

    /* The topology: a control point per residue at its CA, pieces flagged */
    md_temp_scope_t temp = md_temp_begin();
    md_render_mol_plan_t plan;
    ASSERT_TRUE(md_render_mol_plan(&plan, &sys, NULL, temp));
    const md_render_mol_topo_t* L = &plan.layout;
    const uint32_t ncp = L->backbone_count;
    ASSERT_GT(ncp, 50u);
    {
        const uint8_t*  img   = (const uint8_t*)plan.data;
        const uint32_t* order = (const uint32_t*)(img + L->atom_order);
        const uint32_t* idx   = (const uint32_t*)(img + L->backbone_index);
        const uint32_t* seg   = (const uint32_t*)(img + L->backbone_segment);
        const uint32_t* chain = (const uint32_t*)(img + L->backbone_chain);
        int bad = 0;
        for (uint32_t j = 0; j < ncp; ++j) {
            uint32_t ca = order[idx[j] & 0x3FFFFFFFu];
            if (seg[j] >= PB->segment.count || (uint32_t)PB->segment.atoms[seg[j]].ca != ca) bad++;
        }
        EXPECT_EQ(0, bad);
        for (uint32_t c = 0; c < L->backbone_chain_count; ++c) {
            EXPECT_TRUE((idx[chain[c]] & MD_RENDER_SPLINE_FIRST) != 0);
            EXPECT_TRUE((idx[chain[c + 1] - 1] & MD_RENDER_SPLINE_LAST) != 0);
            EXPECT_GE(chain[c + 1] - chain[c], 2u);
        }
        printf("    %u residues in %u pieces, of %u segments\n", ncp, L->backbone_chain_count, (uint32_t)PB->segment.count);
    }
    md_gpu_mem_t mem = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, plan.size);
    md_render_mol_topo_t topo;
    ASSERT_TRUE(md_render_mol_write(fx.gfx, &topo, &plan, mem.gpu, content));
    md_temp_end(temp);

    /* Positions, the update (axes), the secondary structure, the profile */
    const uint32_t na = (uint32_t)sys.atom.count;
    md_gpu_addr_t P = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, (uint64_t)na * 12).gpu;
    ASSERT_TRUE(md_render_mol_upload_state(fx.gfx, P, &state, content));
    md_render_mol_atoms_t atoms = {0};
    atoms.position = P;
    md_render_mol_t mol = md_render_mol_create(fx.dev, NULL);
    ASSERT_TRUE(mol != NULL);
    md_gpu_addr_t dyn = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, topo.dynamic_size).gpu;
    ASSERT_TRUE(md_render_mol_update(mol, fx.gfx, &topo, dyn, &atoms, 0, content));
    ASSERT_TRUE(md_system_state_backbone_compute(&state, &sys));
    const md_secondary_structure_t* ss = md_system_state_secondary_structure(&state, &sys);
    ASSERT_TRUE(ss != NULL);
    md_gpu_addr_t SS = f_buf(&fx, ss, PB->segment.count * sizeof(int32_t));
    md_gpu_addr_t prof = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, md_render_mol_cartoon_size(&topo)).gpu;
    const float scale[3] = { 1.0f, 1.0f, 1.0f };
    ASSERT_TRUE(md_render_mol_cartoon_profile(mol, fx.gfx, &topo, prof, SS, scale, content));

    /* The profile: md_gl's shapes inside runs of one structure, and parts
       meeting without a step (md_render would cap them) */
    {
        float*    pf  = (float*)read_block(&fx, prof, (uint64_t)ncp * 16);
        uint32_t* idx = (uint32_t*)read_block(&fx, topo.backbone_index, (uint64_t)ncp * 4);
        uint32_t* seg = (uint32_t*)read_block(&fx, topo.backbone_segment, (uint64_t)ncp * 4);
        ASSERT_TRUE(pf && idx && seg);
        int runs[3] = { 0, 0, 0 }, bad = 0, steps = 0;
        for (uint32_t j = 0; j < ncp; ++j) {
            if (!(idx[j] & MD_RENDER_SPLINE_LAST) && (pf[4 * j + 1] != pf[4 * (j + 1)] || pf[4 * j + 3] != pf[4 * (j + 1) + 2])) steps++;
            if (j < 2 || j + 2 >= ncp) continue;
            int cls[5];
            bool same = true;
            for (int k = -2; k <= 2; ++k) {
                if (k < 0 && (idx[j + k + 1] & MD_RENDER_SPLINE_FIRST)) same = false;
                if (k > 0 && (idx[j + k - 1] & MD_RENDER_SPLINE_LAST)) same = false;
                int t = ss[seg[j + k]];
                cls[k + 2] = (t >= MD_SECONDARY_STRUCTURE_HELIX_310 && t <= MD_SECONDARY_STRUCTURE_HELIX_PI) ? 1 : (t == MD_SECONDARY_STRUCTURE_BETA_SHEET || t == MD_SECONDARY_STRUCTURE_BETA_BRIDGE) ? 2 : 0;
                if (cls[k + 2] != cls[0]) same = false;
            }
            if (!same) continue;
            const float want[3][2] = { { 0.0f, 0.2f }, { 1.0f, 0.2f }, { 1.4f, 0.1f } };
            const float* w = want[cls[0]];
            runs[cls[0]]++;
            if (fabsf(pf[4 * j] - w[0]) > 1e-6f || fabsf(pf[4 * j + 1] - w[0]) > 1e-6f || fabsf(pf[4 * j + 2] - w[1]) > 1e-6f || fabsf(pf[4 * j + 3] - w[1]) > 1e-6f) {
                if (bad < 4) printf("    control point %u (class %d): %g %g %g %g\n", j, cls[0], pf[4 * j], pf[4 * j + 1], pf[4 * j + 2], pf[4 * j + 3]);
                bad++;
            }
        }
        EXPECT_EQ(0, bad);
        EXPECT_EQ(0, steps);
        EXPECT_GT(runs[0], 5);
        EXPECT_GT(runs[1] + runs[2], 5);
        printf("    profiles inside runs: coil %d, helix %d, sheet %d\n", runs[0], runs[1], runs[2]);
        free(pf); free(idx); free(seg);
    }

    /* The axes: unit, across the trace, and (the first update decides from
       the structure alone) each within 90 degrees of the one before */
    md_render_mol_rep_t rep = {0};
    rep.key  = 7;
    rep.type = MD_RENDER_MOL_REP_CARTOON;
    rep.params.cartoon.profile = prof;
    rep.color.mode = MD_RENDER_COLOR_UNIFORM;
    rep.color.uniform = 0xFF40A0E0u;
    md_render_batch_t batch;
    ASSERT_EQ(1u, md_render_mol_emit(&topo, dyn, &atoms, &rep, &batch, 1));
    EXPECT_EQ(MD_RENDER_SPLINES, batch.kind);
    EXPECT_EQ(ncp, batch.count);
    const md_gpu_addr_t axes = batch.style.splines.axis;
    {
        float* ax = (float*)read_block(&fx, axes, (uint64_t)ncp * 16);
        uint32_t* idx = (uint32_t*)read_block(&fx, topo.backbone_index, (uint64_t)ncp * 4);
        uint32_t* order = (uint32_t*)read_block(&fx, topo.atom_order, (uint64_t)topo.point_count * 4);
        int bad = 0;
        for (uint32_t j = 0; j < ncp; ++j) {
            const float* a = ax + 4 * j;
            float l = sqrtf(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
            uint32_t jp = (idx[j] & MD_RENDER_SPLINE_FIRST) ? j : j - 1, jn = (idx[j] & MD_RENDER_SPLINE_LAST) ? j : j + 1;
            const float* xp = &state.xyz[order[idx[jp] & 0x3FFFFFFFu]].x;
            const float* xn = &state.xyz[order[idx[jn] & 0x3FFFFFFFu]].x;
            float t[3] = { xn[0] - xp[0], xn[1] - xp[1], xn[2] - xp[2] };
            float tl = sqrtf(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]);
            float across = fabsf(a[0] * t[0] + a[1] * t[1] + a[2] * t[2]) / (tl * l);
            float prev = 1.0f;
            if (!(idx[j] & MD_RENDER_SPLINE_FIRST)) prev = a[0] * ax[4 * j - 4] + a[1] * ax[4 * j - 3] + a[2] * ax[4 * j - 2];
            if (fabsf(l - 1.0f) > 1e-4f || across > 1e-3f || prev < 0.0f) bad++;
        }
        EXPECT_EQ(0, bad);
        free(ax); free(idx); free(order);
    }

    /* Drawn: an orthographic view of the whole protein from +z; the curve
       point nearest the viewer (no part in front of it) picks as its residue */
    float lo[3] = { 1e30f, 1e30f, 1e30f }, hi[3] = { -1e30f, -1e30f, -1e30f };
    for (uint32_t i = 0; i < na; ++i) for (int k = 0; k < 3; ++k) { lo[k] = fminf(lo[k], (&state.xyz[i].x)[k]); hi[k] = fmaxf(hi[k], (&state.xyz[i].x)[k]); }
    const uint32_t W = 256, H = 256;
    float ext = fmaxf(hi[0] - lo[0], hi[1] - lo[1]) * 0.55f + 2.0f;
    cam_t c = { .ortho = true, .half_h = ext, .aspect = 1.0f, .n = 0.0f, .far_ = 1000.0f };
    md_render_camera_t cam = make_camera(&c);
    float ctr[3] = { 0.5f * (lo[0] + hi[0]), 0.5f * (lo[1] + hi[1]), hi[2] + 20.0f };
    cam.view[12] = -ctr[0]; cam.view[13] = -ctr[1]; cam.view[14] = -ctr[2];
    md_render_frame_t f = frame_of(W, H, cam, &batch, 1);
    memcpy(f.content, content, sizeof(content));
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    md_gpu_stream_sync(fx.gfx);
    md_render_stats_t st;
    ASSERT_TRUE(md_render_view_stats(fx.view, &st));
    EXPECT_GT(st.hardware[0] + st.hardware[1], ncp / 2);
    {
        uint32_t* idx = (uint32_t*)read_block(&fx, topo.backbone_index, (uint64_t)ncp * 4);
        uint32_t* order = (uint32_t*)read_block(&fx, topo.atom_order, (uint64_t)topo.point_count * 4);
        uint32_t* seg = (uint32_t*)read_block(&fx, topo.backbone_segment, (uint64_t)ncp * 4);
        uint32_t best = 0;
        float bz = -1e30f, bp[3] = { 0, 0, 0 };
        for (uint32_t j = 0; j < ncp; ++j) {
            if (idx[j] & (MD_RENDER_SPLINE_FIRST | MD_RENDER_SPLINE_LAST)) continue;
            const float* a = &state.xyz[order[idx[j - 1] & 0x3FFFFFFFu]].x;
            const float* b = &state.xyz[order[idx[j] & 0x3FFFFFFFu]].x;
            const float* d = &state.xyz[order[idx[j + 1] & 0x3FFFFFFFu]].x;
            float q[3];
            for (int k = 0; k < 3; ++k) q[k] = (a[k] + 4.0f * b[k] + d[k]) / 6.0f;
            if (q[2] > bz) { bz = q[2]; best = j; memcpy(bp, q, sizeof(q)); }
        }
        uint32_t px = (uint32_t)(((bp[0] - ctr[0]) / ext * 0.5f + 0.5f) * W);
        uint32_t py = (uint32_t)((0.5f - (bp[1] - ctr[1]) / ext * 0.5f) * H);
        md_render_pick_t pk;
        md_render_mol_hit_t hit;
        ASSERT_TRUE(f_pick(&fx, px, py, &pk));
        ASSERT_TRUE(md_render_mol_decode_pick(&pk, &hit));
        EXPECT_EQ(MD_RENDER_MOL_HIT_BACKBONE, hit.kind);
        EXPECT_EQ(7u, (unsigned)hit.rep_key);
        EXPECT_EQ(seg[best], hit.segment);
        EXPECT_EQ((uint32_t)PB->segment.atoms[seg[best]].ca, hit.atom);
        EXPECT_EQ(UINT32_MAX, hit.bond);
        free(idx); free(order); free(seg);
    }

    /* Sides over a jittered trajectory (0.15 A per coordinate, 40 steps):
       with memory no axis turns over against its last; decided afresh every
       update (the reset), pairs near 90 degrees apart turn over against each
       other (a ribbon's half turn between two residues, what is seen; a
       whole piece turned over is not, and is not counted) */
    {
        md_float3_t* xyz = malloc((size_t)na * sizeof(md_float3_t));
        uint32_t* idx_host = (uint32_t*)read_block(&fx, topo.backbone_index, (uint64_t)ncp * 4);
        float* last = (float*)read_block(&fx, axes, (uint64_t)ncp * 16);
        int turned[2] = { 0, 0 };
        for (int mode = 0; mode < 2; ++mode) {
            uint32_t rng = 12345u;
            ASSERT_TRUE(md_render_mol_update(mol, fx.gfx, &topo, dyn, &atoms, MD_RENDER_MOL_UPDATE_RESET_BACKBONE, content));
            free(last);
            last = (float*)read_block(&fx, axes, (uint64_t)ncp * 16);
            for (int step = 0; step < 40; ++step) {
                for (uint32_t i = 0; i < na; ++i) {
                    xyz[i].x = state.xyz[i].x + 0.15f * lcg_gauss(&rng);
                    xyz[i].y = state.xyz[i].y + 0.15f * lcg_gauss(&rng);
                    xyz[i].z = state.xyz[i].z + 0.15f * lcg_gauss(&rng);
                }
                md_gpu_upload(fx.gfx, P, xyz, (size_t)na * 12);
                ASSERT_TRUE(md_render_mol_update(mol, fx.gfx, &topo, dyn, &atoms, mode ? MD_RENDER_MOL_UPDATE_RESET_BACKBONE : 0, content));
                float* now = (float*)read_block(&fx, axes, (uint64_t)ncp * 16);
                bool prev_over = false;
                for (uint32_t j = 0; j < ncp; ++j) {
                    float d = now[4 * j] * last[4 * j] + now[4 * j + 1] * last[4 * j + 1] + now[4 * j + 2] * last[4 * j + 2];
                    bool over = d < 0.0f;
                    if (mode == 0 && over) turned[0]++;
                    if (mode == 1 && !(idx_host[j] & MD_RENDER_SPLINE_FIRST) && over != prev_over) turned[1]++;
                    prev_over = over;
                }
                free(last);
                last = now;
            }
        }
        free(last); free(xyz); free(idx_host);
        printf("    in 40 jittered updates: %d axes turned over with memory, %d pairs against each other decided afresh\n", turned[0], turned[1]);
        EXPECT_EQ(0, turned[0]);
        EXPECT_GT(turned[1], 0);
    }

    md_gpu_free(fx.gfx, P);
    md_gpu_free(fx.gfx, dyn);
    md_gpu_free(fx.gfx, prof);
    md_gpu_free(fx.gfx, topo.base);
    md_render_mol_destroy(mol);
    md_system_free(&sys);
    f_close(&fx);
}

#endif
