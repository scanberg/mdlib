#include "utest.h"

#include <md_render.h>
#include <md_render_mol.h>
#include <md_render_debug.h>
#include <md_compose.h>
#include <md_post.h>

#include <md_system.h>
#include <md_system_state.h>
#include <md_pdb.h>
#include <md_gro.h>
#include <md_mmcif.h>
#include <md_infer.h>
#include <md_backbone.h>
#include <md_ring.h>
#include <md_bond.h>
#include <md_chem.h>
#include <md_xyz.h>
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
    md_render_ctx_t  r;
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

/* Marks of n points from bitset words over them (NULL: none set) */
static void f_marks_bytes(uint8_t* m, const uint32_t* s0, const uint32_t* s1, uint32_t n) {
    for (uint32_t i = 0; i < n; ++i) {
        m[i] = (uint8_t)((s0 && ((s0[i >> 5] >> (i & 31)) & 1u) ? MD_RENDER_MARK_STATE0 : 0u) |
                         (s1 && ((s1[i >> 5] >> (i & 31)) & 1u) ? MD_RENDER_MARK_STATE1 : 0u));
    }
}
static void f_marks_upload(fix_t* f, md_gpu_addr_t dst, const uint32_t* s0, const uint32_t* s1, uint32_t n) {
    uint8_t* m = malloc(n);
    f_marks_bytes(m, s0, s1, n);
    md_gpu_upload(f->gfx, dst, m, n);
    free(m);
}
static md_gpu_addr_t f_marks(fix_t* f, const uint32_t* s0, const uint32_t* s1, uint32_t n) {
    md_gpu_addr_t a = md_gpu_malloc(f->gfx, MD_GPU_MEM_DEVICE, n).gpu;
    if (a) f_marks_upload(f, a, s0, s1, n);
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

/* Radii from a table: float4 points whose w is a type, an index into it;
   the spheres as drawn are those of the table's radii, and a type past the
   table's end is radius 0 (not drawn) */
UTEST(render, sphere_radius_table) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 64, H = 64;
    cam_t c = { .ortho = true, .half_h = 3.0f, .aspect = 1.0f, .n = 0.0f, .far_ = 50.0f };
    md_render_camera_t cam = make_camera(&c);
    const double ctr[2][3] = { { 0.5, 0.2, -10 }, { -1.0, -1.0, -20 } };
    const double rad[2]    = { 1.2, 2.5 };
    const float table[3]   = { 0.0f, 2.5f, 1.2f };
    const uint32_t type[3] = { 2u, 1u, 7u };   /* the third: past the table */
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
    md_render_debug_view(hw, MD_RENDER_DEBUG_NO_COMPUTE_RASTER, 0.0f);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, hw, &fh, NULL));
    md_render_debug_view(fx.view, 0, 1000.0f);
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

/* Every pixel of a view of splines against their sweep: empty where the rays
   miss it, its surface where they hit, nothing of it in front, the part
   that is nearest where that is clear */
static void sweep_check_pixels(const cam_t* c, uint32_t W, uint32_t H, const sweep_t* sw, const uint32_t* id, const float* dp,
                               int* bad, int* covered, int* judged_parts, double* worst_surface) {
    *bad = 0; *covered = 0; *judged_parts = 0; *worst_surface = 0;
    for (uint32_t y = 0; y < H; ++y) {
        for (uint32_t x = 0; x < W; ++x) {
            ray_t r = pixel_ray(c, W, H, x + 0.5, y + 0.5);
            uint32_t ix = id[2 * (y * W + x)], iy = id[2 * (y * W + x) + 1];
            double clr = sweep_ray_clearance(sw, r, 100.0);
            /* a pixel's size in world at depth z, from the ray through the scene's middle */
            double px_at15 = 15.0 * 2.0 * c->tan_half / H;
            if (iy == 0) {
                if (clr < -0.75 * px_at15) { (*bad)++; if (*bad < 8) printf("    pixel (%u, %u): empty, the sweep %.3f px inside\n", x, y, -clr / px_at15); }
                continue;
            }
            (*covered)++;
            if (clr > 0.75 * px_at15) { (*bad)++; if (*bad < 8) printf("    pixel (%u, %u): covered, the sweep %.3f px away\n", x, y, clr / px_at15); continue; }
            double zv = -c->n / dp[y * W + x], tp = zv / r.d[2];
            double P[3] = { r.o[0] + r.d[0] * tp, r.o[1] + r.d[1] * tp, r.o[2] + r.d[2] * tp };
            double pxw = -P[2] * 2.0 * c->tan_half / H;
            int who; double gap;
            double sd = sweep_sdf(sw, P, &who, &gap);
            *worst_surface = fmax(*worst_surface, fabs(sd) / pxw);
            if (fabs(sd) > 0.5 * pxw) { (*bad)++; if (*bad < 8) printf("    pixel (%u, %u): %.3f px off the surface\n", x, y, sd / pxw); continue; }
            if (gap > pxw) {
                (*judged_parts)++;
                if (ix != (uint32_t)who) { (*bad)++; if (*bad < 8) printf("    pixel (%u, %u): part %u, expected %d\n", x, y, ix, who); continue; }
            }
            /* nothing of the sweep in front: the ray before P stays outside but for a pixel */
            double tt = 0.0;
            for (int it = 0; it < 400 && tt < tp; ++it) {
                double q[3] = { r.o[0] + r.d[0] * tt, r.o[1] + r.d[1] * tt, r.o[2] + r.d[2] * tt };
                double d = sweep_sdf(sw, q, &who, &gap);
                double qw = -q[2] * 2.0 * c->tan_half / H;
                if (d < -0.5 * qw) { (*bad)++; if (*bad < 8) printf("    pixel (%u, %u): the sweep %.3f px in front\n", x, y, (tp - tt) / qw); break; }
                tt += fmax(fabs(d), 0.25 * qw);
            }
        }
    }
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
    md_render_debug_view(fx.view, MD_RENDER_DEBUG_NO_COMPUTE_RASTER, 0.0f);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));

    sweep_t sw;
    sweep_build(&sw, (const double (*)[3])vpos, (const double (*)[3])vaxis, profile, index, visible, N);
    md_render_targets_t tg = md_render_view_targets(fx.view);
    uint32_t* id = malloc(W * H * 8);
    float*    dp = malloc(W * H * 4);
    ASSERT_TRUE(f_read(&fx, tg.id, id, W * H * 8) && f_read(&fx, tg.depth, dp, W * H * 4));
    int bad = 0, covered = 0, judged_parts = 0;
    double worst_surface = 0;    /* in pixels */
    sweep_check_pixels(&c, W, H, &sw, id, dp, &bad, &covered, &judged_parts, &worst_surface);
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
    md_render_debug_view(cr, 0, 1000.0f);
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

    /* Through its control points: each moved first, as md_render does
       (c = (8 q - q- - q+) / 6, the ends staying), and the sweep of those */
    {
        double tpos[N][3], before = 0, after = 0;
        for (int k = 0; k < N; ++k) {
            bool end = (index[k] & (MD_RENDER_SPLINE_FIRST | MD_RENDER_SPLINE_LAST)) != 0;
            for (int i = 0; i < 3; ++i) tpos[k][i] = end ? vpos[k][i] : (8.0 * vpos[k][i] - vpos[k - 1][i] - vpos[k + 1][i]) / 6.0;
        }
        for (int k = 1; k < NA - 1; ++k) {       /* the helix: how far its knots are from its points */
            double d0[3], d1[3];
            for (int i = 0; i < 3; ++i) {
                d0[i] = (vpos[k - 1][i] + 4.0 * vpos[k][i] + vpos[k + 1][i]) / 6.0 - vpos[k][i];
                d1[i] = (tpos[k - 1][i] + 4.0 * tpos[k][i] + tpos[k + 1][i]) / 6.0 - vpos[k][i];
            }
            before = fmax(before, sqrt(dot3(d0, d0)));
            after  = fmax(after, sqrt(dot3(d1, d1)));
        }
        printf("    through: the helix's knots off its points by %.3f A at most, %.3f near them\n", after, before);
        EXPECT_LT(after, 0.45 * before);
        b.style.splines.through = true;
        f.content[MD_RENDER_CONTENT_GEOMETRY]++;
        ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
        sweep_t tsw;
        sweep_build(&tsw, (const double (*)[3])tpos, (const double (*)[3])vaxis, profile, index, visible, N);
        uint32_t* tid = malloc(W * H * 8);
        float*    tdp = malloc(W * H * 4);
        ASSERT_TRUE(f_read(&fx, tg.id, tid, W * H * 8) && f_read(&fx, tg.depth, tdp, W * H * 4));
        int tbad = 0, tcov = 0, tjudged = 0;
        double tworst = 0;
        sweep_check_pixels(&c, W, H, &tsw, tid, tdp, &tbad, &tcov, &tjudged, &tworst);
        printf("    through: %d pixels covered, %d with a clear part; the farthest %.3f px from the sweep\n", tcov, tjudged, tworst);
        EXPECT_EQ(0, tbad);
        EXPECT_GT(tcov, covered * 3 / 4);
        EXPECT_GT(tjudged, tcov / 2);
        free(tid); free(tdp); free(tsw.part);
        b.style.splines.through = false;
        f.content[MD_RENDER_CONTENT_GEOMETRY]++;
    }

    /* every other part of the helix selected: its mask where the raster
       shows them, its outline around more (where the turn in front covers
       the turn behind) */
    uint32_t sel = 0x555u & ((1u << NA) - 1u);
    b.points.marks = f_marks(&fx, &sel, NULL, 32);
    f.content[MD_RENDER_CONTENT_STATE]++;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    /* against this frame's ids, which the mask is made from: it rasterised
       again (the through pass changed the geometry and back), and where two
       parts meet their surfaces tie, so which one a seam pixel shows follows
       the order they were drawn in, which the cull's appends do not fix (the
       HD 530 showed part 7 where the first frame had its neighbour) */
    ASSERT_TRUE(f_read(&fx, tg.id, id, W * H * 8));
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
    b.points.marks = 0;
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
        /* by control points, in front and behind; each its own (k * 7, the
           first take, put k and k + 12 in one place: N is 28, and which of
           two equal depths wins is the draw order's, no two views' alike) */
        int e = (k * 5) % N;
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

    /* every spline pixel coloured, by the heavy resolve over the tiles the
       main one listed: a target of whole tiles, and one whose last row and
       column of tiles are partial */
    {
        md_render_batch_t bc = b;
        bc.color.mode    = MD_RENDER_COLOR_UNIFORM;
        bc.color.uniform = 0xFF4080C0u;
        md_render_batch_t both_c[2] = { bs, bc };
        const uint32_t sz[2][2] = { { W, H }, { W - 3, H - 5 } };
        for (int k = 0; k < 2; ++k) {
            const uint32_t w = sz[k][0], h = sz[k][1];
            md_render_view_t vc = md_render_view_create(fx.r, "colours");
            md_render_frame_t fc2 = frame_of(w, h, cam, both_c, 2);
            ASSERT_TRUE(md_render_view_draw(fx.gfx, vc, &fc2, NULL));
            md_gpu_stream_sync(fx.gfx);
            uint32_t* idc = malloc((size_t)w * h * 8);
            uint32_t* cc  = malloc((size_t)w * h * 4);
            md_render_targets_t t = md_render_view_targets(vc);
            ASSERT_TRUE(f_read(&fx, t.id, idc, (size_t)w * h * 8) && f_read(&fx, t.color, cc, (size_t)w * h * 4));
            int spl = 0, sph = 0, wrong = 0;
            for (uint32_t i = 0; i < w * h; ++i) {
                uint32_t slot = idc[2 * i + 1] >> 20;
                if (slot == 0) continue;
                uint32_t want = slot == 2 ? 0xFF4080C0u : 0xFF30C030u;
                if (slot == 2) spl++; else sph++;
                if (!near_rgba(cc[i], want)) {
                    if (wrong < 4) printf("    %ux%u, pixel (%u, %u) of %s: colour %08x\n", w, h, i % w, i / w, slot == 2 ? "the splines" : "a sphere", cc[i]);
                    wrong++;
                }
            }
            EXPECT_EQ(0, wrong);
            EXPECT_GT(spl, 1000);
            EXPECT_GT(sph, 200);
            printf("    %ux%u: %d spline pixels and %d sphere pixels coloured\n", w, h, spl, sph);
            free(idc); free(cc);
            md_render_view_destroy(fx.gfx, vc);
        }
    }
    free(id); free(dp); free(sw.part);
    f_close(&fx);
}

/* ---- Polygons ------------------------------------------------------------------ */

/* A polygon as md_render.h defines it, in double: the corners' centroid,
   Newell's normal, the corners projected onto that plane; rounded by r. On a
   line: its edges alone. */
typedef struct poly_ref_t {
    double q[8][3], c[3], n[3], r, bound;
    double h;                            /* the plate's half thickness: r flat, less a dish */
    int    corners;
    bool   planar;
    uint32_t key, linear;                /* its batch, and its id's element (instance * count + element) */
    uint32_t point[8];                   /* its corners' points: their attributes */
} poly_ref_t;

static void cross3(const double* a, const double* b, double* o) {
    o[0] = a[1] * b[2] - a[2] * b[1]; o[1] = a[2] * b[0] - a[0] * b[2]; o[2] = a[0] * b[1] - a[1] * b[0];
}

static void poly_ref_build(poly_ref_t* p, const double (*q)[3], int n, double r) {
    p->corners = n;
    p->r = r;
    p->h = r;
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
    p->planar = nl > 1e-6 * sz * sz;
    for (int i = 0; i < 3; ++i) p->n[i] = p->planar ? nv[i] / nl : 0.0;
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
    if (!p->planar) return false;
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

/* The signed distance to the polygon's solid. Flat: to the polygon, less
   the radius. A dish (h < r): the nearer of its rims (the capsules of its
   edges) and its plate (the prism of the polygon, h either side), and how
   near those two are (crease; 1e30 when flat) */
static double poly_ref_sdf(const poly_ref_t* p, const double x[3], double* crease) {
    double cp[3];
    if (crease) *crease = 1e30;
    if (p->h >= p->r || !p->planar) return poly_ref_dist(p, x, cp) - p->r;
    double rim = 1e30, edge = 1e30;
    double d[3] = { x[0] - p->c[0], x[1] - p->c[1], x[2] - p->c[2] }, z = dot3(d, p->n);
    double xp[3] = { x[0] - p->n[0] * z, x[1] - p->n[1] * z, x[2] - p->n[2] * z };
    for (int k = 0; k < p->corners; ++k) {
        rim  = fmin(rim, pt_seg_dist(x, p->q[k], p->q[(k + 1) % p->corners]) - p->r);
        edge = fmin(edge, pt_seg_dist(xp, p->q[k], p->q[(k + 1) % p->corners]));
    }
    double wx = poly_ref_inside(p, xp) ? -edge : edge, wy = fabs(z) - p->h;
    double plate = fmin(fmax(wx, wy), 0.0) + sqrt(fmax(wx, 0.0) * fmax(wx, 0.0) + fmax(wy, 0.0) * fmax(wy, 0.0));
    if (crease) *crease = fabs(rim - plate);
    return fmin(rim, plate);
}

/* The signed distance to the nearest polygon, which, and how much farther the next one is */
static double polys_sdf(const poly_ref_t* P, int count, const double x[3], int* who, double* gap) {
    double best = 1e30, second = 1e30;
    *who = -1;
    for (int k = 0; k < count; ++k) {
        double d = poly_ref_sdf(&P[k], x, NULL);
        if (d < best) { second = best; best = d; *who = k; }
        else if (d < second) second = d;
    }
    *gap = second - best;
    return best;
}

/* How far the ray (t in [0, tmax]) passes outside every polygon: < 0, through
   one. A dish: its rims, or how far inside the polygon it crosses the plate's
   middle, at most h (a ray into the plate either crosses that or leaves
   through a rim) */
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
        if (p->planar && p->h < p->r) {
            double rim = 1e30, plate = 1e30;
            for (int j = 0; j < p->corners; ++j) rim = fmin(rim, seg_seg_dist(r.o, q, p->q[j], p->q[(j + 1) % p->corners]) - p->r);
            if ((h0 > 0) != (h1 > 0)) {
                double s = h0 / (h0 - h1), x[3] = { r.o[0] + (q[0] - r.o[0]) * s, r.o[1] + (q[1] - r.o[1]) * s, r.o[2] + (q[2] - r.o[2]) * s };
                double edge = 1e30;
                for (int j = 0; j < p->corners; ++j) edge = fmin(edge, pt_seg_dist(x, p->q[j], p->q[(j + 1) % p->corners]));
                plate = poly_ref_inside(p, x) ? -fmin(edge, p->h) : edge;
            }
            best = fmin(best, fmin(rim, plate));
            continue;
        }
        if (p->planar && (h0 > 0) != (h1 > 0)) {
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

/* Every pixel of a view of polygons against their solids: empty where the
   rays miss them, covered where they hit (flat ones: a dish's coverage is
   its surface's), on the surface, nothing in front; where one polygon is
   clearly the nearest its id and its normal */
static void poly_check_pixels(const cam_t* c, uint32_t W, uint32_t H, const poly_ref_t* ref, int nref, const uint32_t keys[4], bool dished,
                              const uint32_t* id, const float* dp, const uint32_t* nb,
                              int* bad, int* covered, int* judged, double* worst, double* errs, int* nerr, int* per_ref) {
    *bad = 0; *covered = 0; *judged = 0; *worst = 0; *nerr = 0;
    for (int k = 0; k < nref; ++k) per_ref[k] = 0;
    for (uint32_t y = 0; y < H; ++y) {
        for (uint32_t x = 0; x < W; ++x) {
            uint32_t i = y * W + x;
            ray_t r = pixel_ray(c, W, H, x + 0.5, y + 0.5);
            uint32_t ix = id[2 * i], iy = id[2 * i + 1];
            double clr = polys_ray_clearance(ref, nref, r, 100.0);
            double px15 = 15.0 * 2.0 * c->tan_half / H;
            if (iy == 0) {
                if (clr < -0.5 * px15) { (*bad)++; if (*bad < 8) printf("    pixel (%u, %u): empty, a polygon %.3f px inside\n", x, y, -clr / px15); }
                continue;
            }
            (*covered)++;
            if (!dished && clr > 0.5 * px15) { (*bad)++; if (*bad < 8) printf("    pixel (%u, %u): covered, the polygons %.3f px away\n", x, y, clr / px15); continue; }
            double zv = -c->n / dp[i], tp = zv / r.d[2];
            double P[3] = { r.o[0] + r.d[0] * tp, r.o[1] + r.d[1] * tp, r.o[2] + r.d[2] * tp };
            double pxw = -P[2] * 2.0 * c->tan_half / H;
            int who; double gap;
            double sd = polys_sdf(ref, nref, P, &who, &gap);
            *worst = fmax(*worst, fabs(sd) / pxw);
            if (fabs(sd) > 0.05 * pxw) { (*bad)++; if (*bad < 8) printf("    pixel (%u, %u): %.3f px off the surface\n", x, y, sd / pxw); continue; }
            if (gap > pxw) {
                (*judged)++;
                per_ref[who]++;
                uint32_t slot = (iy >> 20) - 1u;
                if (slot >= 4 || keys[slot] != ref[who].key || ix != ref[who].linear) {
                    (*bad)++;
                    if (*bad < 8) printf("    pixel (%u, %u): batch slot %u element %u, expected key %u element %u\n", x, y, slot, ix, ref[who].key, ref[who].linear);
                    continue;
                }
                /* the normal: the distance's gradient, away from a dish's crease */
                double nr[3], ng[3], crease;
                poly_ref_sdf(&ref[who], P, &crease);
                for (int k = 0; k < 3; ++k) {
                    double a[3] = { P[0], P[1], P[2] }, b[3] = { P[0], P[1], P[2] };
                    a[k] += 1e-6; b[k] -= 1e-6;
                    nr[k] = (poly_ref_sdf(&ref[who], a, NULL) - poly_ref_sdf(&ref[who], b, NULL)) / 2e-6;
                }
                double nl = sqrt(dot3(nr, nr));
                if (crease > 0.5 * pxw && fabs(dot3(nr, r.d)) > 0.3 * nl) {
                    for (int k = 0; k < 3; ++k) nr[k] /= nl;
                    oct_decode(nb[i], ng);
                    double cs = dot3(nr, ng);
                    errs[(*nerr)++] = acos(cs > 1 ? 1 : cs < -1 ? -1 : cs) * 180.0 / 3.14159265358979;
                }
            }
            /* nothing in front of it */
            double tt = 0.0;
            for (int it = 0; it < 400 && tt < tp; ++it) {
                double qv[3] = { r.o[0] + r.d[0] * tt, r.o[1] + r.d[1] * tt, r.o[2] + r.d[2] * tt };
                double d = polys_sdf(ref, nref, qv, &who, &gap);
                double qw = -qv[2] * 2.0 * c->tan_half / H;
                if (d < -0.05 * qw) { (*bad)++; if (*bad < 8) printf("    pixel (%u, %u): a polygon %.3f px in front\n", x, y, (tp - tt) / qw); break; }
                tt += fmax(fabs(d), 0.05 * qw);
            }
        }
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
    EXPECT_FALSE(ref[NH - 1 + 2 * NP + 1].planar);        /* the triangle on a line */

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
    poly_check_pixels(&c, W, H, ref, nref, keys, false, id, dp, nb, &bad, &covered, &judged, &worst, errs, &nerr, per_ref);
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
    bp[0].points.marks = f_marks(&fx, hsel, NULL, (uint32_t)sizeof(hsel) * 8);
    bp[1].points.marks = f_marks(&fx, &psel, NULL, 32);
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
    bp[0].points.marks = 0;
    bp[1].points.marks = 0;

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

    /* dished: the plate half as thick as the rims (md_gl's slabs.plate) */
    {
        for (int k = 0; k < 4; ++k) bp[k].style.polygons.plate = 0.5f;
        for (int k = 0; k < nref; ++k) ref[k].h = 0.5 * ref[k].r;
        f.content[MD_RENDER_CONTENT_GEOMETRY]++;
        ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
        ASSERT_TRUE(f_read(&fx, tg.id, id, W * H * 8) && f_read(&fx, tg.depth, dp, W * H * 4) && f_read(&fx, tg.normal, nb, W * H * 4));
        double* derrs = malloc(W * H * sizeof(double));
        int dbad, dcov, djudged, dnerr, dper[NREF];
        double dworst;
        poly_check_pixels(&c, W, H, ref, nref, keys, true, id, dp, nb, &dbad, &dcov, &djudged, &dworst, derrs, &dnerr, dper);
        qsort(derrs, (size_t)dnerr, sizeof(double), cmp_double);
        double d50 = dnerr ? derrs[dnerr / 2] : 0, d99 = dnerr ? derrs[(dnerr * 99) / 100] : 0;
        printf("    dished: %d pixels covered, %d with a clear polygon; the farthest %.4f px from the solid; normals over %d: median %.3f, 99th %.3f degrees\n",
               dcov, djudged, dworst, dnerr, d50, d99);
        EXPECT_EQ(0, dbad);
        EXPECT_GT(dcov, covered * 9 / 10);
        EXPECT_GT(djudged, dcov * 3 / 4);
        for (int k = 0; k < nref; ++k) EXPECT_GT(dper[k], 20);
        EXPECT_LT(d50, 0.2);
        EXPECT_LT(d99, 1.0);
        free(derrs);
        for (int k = 0; k < 4; ++k) bp[k].style.polygons.plate = 0.0f;
        for (int k = 0; k < nref; ++k) ref[k].h = ref[k].r;
        f.content[MD_RENDER_CONTENT_GEOMETRY]++;
    }

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

/* ---- Dashed capsules ----------------------------------------------------------- */

/* A dashed capsule as md_render.h defines it, in double: the segment between
   where its rim meets the spheres at its ends, dashes of the batch's
   lengths, as many as come closest (at least two), stretched alike to fill it */
typedef struct dash_ref_t {
    double a[3], u[3], len, r;      /* the trimmed segment: start, direction, length; radius */
    double dash, period;
    int    n;                       /* dashes; 0: none (the spheres meet) */
    uint32_t key, linear;
} dash_ref_t;

static void dash_ref_build(dash_ref_t* d, const double pa[3], const double pb[3], double ea, double eb, double r, double D, double G) {
    double ab[3] = { pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2] }, l = sqrt(dot3(ab, ab));
    /* from where the rim meets the sphere of radius ea, eb */
    double ia = sqrt(fmax(ea * ea - r * r, 0.0)), ib = sqrt(fmax(eb * eb - r * r, 0.0));
    d->r = r;
    d->n = 0;
    if (l - ia - ib <= 0.0 || (ia + ib > 0.0 && l - ia - ib <= 0.5 * D)) return;
    for (int i = 0; i < 3; ++i) { d->u[i] = ab[i] / l; d->a[i] = pa[i] + d->u[i] * ia; }
    d->len = l - ia - ib;
    double n = fmax(2.0, floor((d->len + G) / (D + G) + 0.5));
    double k = d->len / (n * D + (n - 1) * G);
    d->n = (int)n; d->dash = D * k; d->period = (D + G) * k;
}

/* The distance to the dashes: the nearest flat-ended cylinder (dash k from k period, dash long) */
static double dash_ref_sdf(const dash_ref_t* d, const double x[3]) {
    if (!d->n) return 1e30;
    double ox[3] = { x[0] - d->a[0], x[1] - d->a[1], x[2] - d->a[2] }, y = dot3(ox, d->u);
    double rp[3] = { ox[0] - d->u[0] * y, ox[1] - d->u[1] * y, ox[2] - d->u[2] * y }, rad = sqrt(dot3(rp, rp));
    double best = 1e30;
    int k0 = (int)floor(y / d->period);
    for (int k = k0 - 1; k <= k0 + 1; ++k) {
        if (k < 0 || k >= d->n) continue;
        double c = k * d->period + 0.5 * d->dash;
        double qx = rad - d->r, qy = fabs(y - c) - 0.5 * d->dash;
        double s = sqrt(fmax(qx, 0) * fmax(qx, 0) + fmax(qy, 0) * fmax(qy, 0)) + fmin(fmax(qx, qy), 0.0);
        best = fmin(best, s);
    }
    if (best >= 1e29) {                    /* beyond the ends: the end dash */
        int k = y < 0 ? 0 : d->n - 1;
        double c = k * d->period + 0.5 * d->dash;
        double qx = rad - d->r, qy = fabs(y - c) - 0.5 * d->dash;
        best = sqrt(fmax(qx, 0) * fmax(qx, 0) + fmax(qy, 0) * fmax(qy, 0)) + fmin(fmax(qx, qy), 0.0);
    }
    return best;
}

/* Whether the ray hits the dashes grown (e > 0) or shrunk (e < 0) by e: each
   dash a finite cylinder, the ray's interval inside the infinite one and
   inside its slab intersected (independent of md_render's way: the next end) */
static bool dash_ref_hits(const dash_ref_t* d, ray_t ry, double e) {
    if (!d->n) return false;
    double oa[3] = { ry.o[0] - d->a[0], ry.o[1] - d->a[1], ry.o[2] - d->a[2] };
    double du = dot3(ry.d, d->u), ou = dot3(oa, d->u);
    double dp[3] = { ry.d[0] - d->u[0] * du, ry.d[1] - d->u[1] * du, ry.d[2] - d->u[2] * du };
    double op[3] = { oa[0] - d->u[0] * ou, oa[1] - d->u[1] * ou, oa[2] - d->u[2] * ou };
    double r = d->r + e, A = dot3(dp, dp), B = dot3(dp, op), C = dot3(op, op) - r * r;
    double t0 = -1e30, t1 = 1e30;
    if (A > 1e-15) {
        double h = B * B - A * C;
        if (h < 0) return false;
        t0 = (-B - sqrt(h)) / A; t1 = (-B + sqrt(h)) / A;
    } else if (C > 0) return false;
    for (int k = 0; k < d->n; ++k) {
        double y0 = k * d->period - e, y1 = k * d->period + d->dash + e, s0 = -1e30, s1 = 1e30;
        if (fabs(du) > 1e-15) {
            s0 = (y0 - ou) / du; s1 = (y1 - ou) / du;
            if (s0 > s1) { double t = s0; s0 = s1; s1 = t; }
        } else if (ou < y0 || ou > y1) continue;
        double lo = fmax(fmax(t0, s0), 0.0), hi = fmin(t1, s1);
        if (lo <= hi) return true;
    }
    return false;
}

UTEST(render, capsules_dashed_match_their_cylinders) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 192, H = 144;
    cam_t c = { .tan_half = 0.4f, .aspect = (float)W / H, .n = 0.1f };
    md_render_camera_t cam = make_camera(&c);
    /* Sticks between points of their own radii: across the view, oblique,
       nearly along it (the ends of the dashes face the camera), short ones,
       one whose spheres meet; a second instance turned and scaled 1.3 */
    enum { NC = 9 };
    const float pos[2 * NC][3] = {
        { -5.0f,  3.0f, -14.0f }, { -1.0f,  3.2f, -14.0f },
        {  0.5f,  3.0f, -13.0f }, {  4.5f,  1.5f, -16.0f },
        { -5.0f,  0.5f, -12.0f }, { -4.6f,  0.9f, -16.5f },     /* nearly along the view */
        { -2.5f,  0.0f, -14.0f }, { -1.6f, -0.2f, -14.0f },     /* short: two dashes */
        {  0.0f,  0.0f, -14.0f }, {  0.6f,  0.0f, -14.0f },     /* the spheres meet */
        {  1.5f, -0.5f, -15.0f }, {  5.5f,  1.0f, -12.0f },
        { -5.5f, -2.5f, -14.0f }, { -1.5f, -3.5f, -15.0f },
        {  0.0f, -2.0f, -14.0f }, {  0.4f, -2.4f, -18.0f },     /* steeply away */
        {  2.0f, -3.0f, -13.5f }, {  5.5f, -2.0f, -14.5f },
    };
    float rad[2 * NC];
    for (int k = 0; k < 2 * NC; ++k) rad[k] = 0.5f + 0.25f * (float)(k % 3);
    const float R = 0.18f, D = 0.35f, G = 0.25f, ES = 0.4f, EA = 0.05f;
    float xf[32] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };
    {   /* about (0, 0, -14): turned 30 degrees about y, scaled 1.3, moved back and down */
        const double ang = 30 * 3.14159265358979 / 180, sc = 1.3, ctr[3] = { 0, 0, -14 }, mv[3] = { 0.5, -1.0, -6.0 };
        float m[9] = { (float)(sc * cos(ang)), 0, (float)(-sc * sin(ang)), 0, (float)sc, 0, (float)(sc * sin(ang)), 0, (float)(sc * cos(ang)) };
        for (int i = 0; i < 9; ++i) xf[16 + 4 * (i / 3) + i % 3] = m[i];
        for (int i = 0; i < 3; ++i) xf[28 + i] = (float)(ctr[i] + mv[i] - (m[i] * ctr[0] + m[3 + i] * ctr[1] + m[6 + i] * ctr[2]));
    }
    md_render_batch_t b = {0};
    b.key = 41; b.kind = MD_RENDER_CAPSULES; b.count = NC; b.points.count = 2 * NC;
    b.points.position = f_buf(&fx, pos, sizeof(pos));
    b.points.radius   = f_buf(&fx, rad, sizeof(rad));
    b.style.capsules.radius = R;
    b.style.capsules.dash_length = D;
    b.style.capsules.gap_length  = G;
    b.style.capsules.dash_ends.radius_scale = ES;
    b.style.capsules.dash_ends.radius_add   = EA;
    b.transform = xf; b.transform_count = 2;
    md_render_frame_t f = frame_of(W, H, cam, &b, 1);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));

    dash_ref_t ref[2 * NC];
    int nref = 0, none = 0;
    for (int in = 0; in < 2; ++in) {
        const float* m = xf + 16 * in;
        const double s = sqrt((double)m[0] * m[0] + (double)m[1] * m[1] + (double)m[2] * m[2]);
        for (int k = 0; k < NC; ++k) {
            double p[2][3];
            for (int e = 0; e < 2; ++e) for (int i = 0; i < 3; ++i) {
                const float* q = pos[2 * k + e];
                p[e][i] = (double)m[i] * q[0] + (double)m[4 + i] * q[1] + (double)m[8 + i] * q[2] + m[12 + i];
            }
            dash_ref_build(&ref[nref], p[0], p[1], (rad[2 * k] * ES + EA) * s, (rad[2 * k + 1] * ES + EA) * s, R * s, D * s, G * s);
            ref[nref].key = 41; ref[nref].linear = (uint32_t)(in * NC + k);
            none += ref[nref].n == 0;
            nref++;
        }
    }
    EXPECT_EQ(2, none);                    /* the one whose spheres meet, twice */
    EXPECT_EQ(2, ref[3].n);                /* the short one: two */

    md_render_targets_t tg = md_render_view_targets(fx.view);
    uint32_t* id = malloc(W * H * 8);
    float*    dp = malloc(W * H * 4);
    uint32_t* nb = malloc(W * H * 4);
    ASSERT_TRUE(f_read(&fx, tg.id, id, W * H * 8) && f_read(&fx, tg.depth, dp, W * H * 4) && f_read(&fx, tg.normal, nb, W * H * 4));
    int bad = 0, covered = 0, judged = 0, nerr = 0, caps = 0;
    double worst = 0;
    double* errs = malloc(W * H * sizeof(double));
    for (uint32_t y = 0; y < H; ++y) {
        for (uint32_t x = 0; x < W; ++x) {
            uint32_t i = y * W + x;
            ray_t r = pixel_ray(&c, W, H, x + 0.5, y + 0.5);
            double px = 14.0 * 2.0 * c.tan_half / H;       /* a pixel at the depth of most */
            bool hit_in = false, hit_out = false;
            for (int k = 0; k < nref; ++k) {
                hit_in  |= dash_ref_hits(&ref[k], r, -0.5 * px);
                hit_out |= dash_ref_hits(&ref[k], r, 0.5 * px);
            }
            uint32_t iy = id[2 * i + 1];
            if (iy == 0) {
                if (hit_in) { bad++; if (bad < 8) printf("    pixel (%u, %u): empty, a dash there\n", x, y); }
                continue;
            }
            covered++;
            if (!hit_out) { bad++; if (bad < 8) printf("    pixel (%u, %u): covered, no dash near\n", x, y); continue; }
            double zv = -c.n / dp[i], tp = zv / r.d[2];
            double P[3] = { r.o[0] + r.d[0] * tp, r.o[1] + r.d[1] * tp, r.o[2] + r.d[2] * tp };
            double pxw = -P[2] * 2.0 * c.tan_half / H;
            int who = -1;
            double best = 1e30, second = 1e30;
            for (int k = 0; k < nref; ++k) {
                double sd = dash_ref_sdf(&ref[k], P);
                if (sd < best) { second = best; best = sd; who = k; } else if (sd < second) second = sd;
            }
            worst = fmax(worst, fabs(best) / pxw);
            if (fabs(best) > 0.05 * pxw) { bad++; if (bad < 8) printf("    pixel (%u, %u): %.3f px off the dashes\n", x, y, best / pxw); continue; }
            if (second - best > pxw) {
                judged++;
                if (id[2 * i] != ref[who].linear) { bad++; if (bad < 8) printf("    pixel (%u, %u): element %u, expected %u\n", x, y, id[2 * i], ref[who].linear); continue; }
                /* the normal: the distance's gradient, away from a dash's rim */
                double g[3], ng[3];
                for (int k = 0; k < 3; ++k) {
                    double a[3] = { P[0], P[1], P[2] }, bq[3] = { P[0], P[1], P[2] };
                    a[k] += 1e-6; bq[k] -= 1e-6;
                    g[k] = (dash_ref_sdf(&ref[who], a) - dash_ref_sdf(&ref[who], bq)) / 2e-6;
                }
                double gl = sqrt(dot3(g, g));
                double ox[3] = { P[0] - ref[who].a[0], P[1] - ref[who].a[1], P[2] - ref[who].a[2] }, yy = dot3(ox, ref[who].u);
                double rp[3] = { ox[0] - ref[who].u[0] * yy, ox[1] - ref[who].u[1] * yy, ox[2] - ref[who].u[2] * yy };
                double radial = sqrt(dot3(rp, rp)), kk = floor(yy / ref[who].period), ya = yy - kk * ref[who].period;
                double edge = fmax(fabs(radial - ref[who].r), fmin(fabs(ya), fabs(ya - ref[who].dash)));   /* the rim: both small */
                if (edge > 0.7 * pxw && fabs(dot3(g, r.d)) > 0.3 * gl) {
                    oct_decode(nb[i], ng);
                    double cs = dot3(g, ng) / gl;
                    errs[nerr++] = acos(cs > 1 ? 1 : cs < -1 ? -1 : cs) * 180.0 / 3.14159265358979;
                    caps += fabs(dot3(g, ref[who].u)) > 0.9 * gl;
                }
            }
        }
    }
    qsort(errs, (size_t)nerr, sizeof(double), cmp_double);
    double p50 = nerr ? errs[nerr / 2] : 0, p99 = nerr ? errs[(nerr * 99) / 100] : 0;
    printf("    %d pixels covered, %d with a clear dash; the farthest %.4f px from the dashes; normals over %d (%d on ends): median %.3f, 99th %.3f degrees\n",
           covered, judged, worst, nerr, caps, p50, p99);
    EXPECT_EQ(0, bad);
    EXPECT_GT(covered, 1000);
    EXPECT_GT(judged, covered * 3 / 4);
    EXPECT_GT(caps, 30);
    EXPECT_LT(p50, 0.2);
    EXPECT_LT(p99, 1.0);
    free(errs);

    /* a pick on the middle of the first dash of the first: that capsule, the nearer end's attribute */
    {
        const dash_ref_t* d = &ref[0];
        double q[3];
        for (int i = 0; i < 3; ++i) q[i] = d->a[i] + d->u[i] * 0.5 * d->dash;
        uint32_t qx = (uint32_t)((q[0] / -q[2] / (c.tan_half * c.aspect) * 0.5 + 0.5) * W);
        uint32_t qy = (uint32_t)((0.5 - q[1] / -q[2] / c.tan_half * 0.5) * H);
        md_render_pick_t pk;
        ASSERT_TRUE(f_pick(&fx, qx, qy, &pk));
        EXPECT_EQ(41u, (unsigned)pk.key);
        EXPECT_EQ(0u, pk.element);
        EXPECT_EQ(0u, pk.attr);
    }
    free(id); free(dp); free(nb);
    f_close(&fx);
}

/* ---- Shaped capsules ------------------------------------------------------- */

/* A capped cone (flat ends) or a round cone (the hull of its end balls),
   radius ra at a to rb at b: the exact signed distance (Quilez's) */
typedef struct {
    double   a[3], b[3], ra, rb;
    bool     flat;
    uint32_t key, linear;
} cone_ref_t;

static double cone_ref_sdf(const cone_ref_t* k, const double p[3]) {
    double ba[3] = { k->b[0] - k->a[0], k->b[1] - k->a[1], k->b[2] - k->a[2] };
    double pa[3] = { p[0] - k->a[0], p[1] - k->a[1], p[2] - k->a[2] };
    double baba = dot3(ba, ba), papa = dot3(pa, pa), y = dot3(pa, ba);
    if (k->flat) {
        double rba = k->rb - k->ra, paba = y / baba;
        double x2 = papa - paba * paba * baba, x = sqrt(x2 > 0 ? x2 : 0);
        double cax = fmax(0.0, x - (paba < 0.5 ? k->ra : k->rb)), cay = fabs(paba - 0.5) - 0.5;
        double kk = rba * rba + baba;
        double f = (rba * (x - k->ra) + paba * baba) / kk;
        f = f < 0 ? 0 : f > 1 ? 1 : f;
        double cbx = x - k->ra - f * rba, cby = paba - f;
        double sg = (cbx < 0.0 && cay < 0.0) ? -1.0 : 1.0;
        return sg * sqrt(fmin(cax * cax + cay * cay * baba, cbx * cbx + cby * cby * baba));
    }
    double rr = k->ra - k->rb, a2 = baba - rr * rr, il2 = 1.0 / baba, z = y - baba;
    double q[3] = { pa[0] * baba - ba[0] * y, pa[1] * baba - ba[1] * y, pa[2] * baba - ba[2] * y };
    double x2 = dot3(q, q), y2 = y * y * baba, z2 = z * z * baba;
    double kk = (rr > 0 ? 1 : rr < 0 ? -1 : 0) * rr * rr * x2;
    if ((z > 0 ? 1 : z < 0 ? -1 : 0) * a2 * z2 > kk) return sqrt(x2 + z2) * il2 - k->rb;
    if ((y > 0 ? 1 : y < 0 ? -1 : 0) * a2 * y2 < kk) return sqrt(x2 + y2) * il2 - k->ra;
    return (sqrt(x2 * a2 * il2) + y * rr) * il2 - k->ra;
}

/* The least distance along the ray (convex: unimodal in t), in pixels at
   its depth (px: a pixel's size per unit depth) */
static double cone_ref_ray_min(const cone_ref_t* k, ray_t r, double px) {
    double c[3] = { 0.5 * (k->a[0] + k->b[0]), 0.5 * (k->a[1] + k->b[1]), 0.5 * (k->a[2] + k->b[2]) };
    double ab[3] = { k->b[0] - k->a[0], k->b[1] - k->a[1], k->b[2] - k->a[2] };
    double R = 0.5 * sqrt(dot3(ab, ab)) + fmax(k->ra, k->rb);
    double oc[3] = { c[0] - r.o[0], c[1] - r.o[1], c[2] - r.o[2] }, tc = dot3(oc, r.d);
    double lo = fmax(tc - R, 0.0), hi = fmax(tc + R, 0.0);
    const double g = 0.5 * (sqrt(5.0) - 1.0);
    for (int it = 0; it < 90; ++it) {
        double t1 = hi - g * (hi - lo), t2 = lo + g * (hi - lo);
        double p1[3] = { r.o[0] + r.d[0] * t1, r.o[1] + r.d[1] * t1, r.o[2] + r.d[2] * t1 };
        double p2[3] = { r.o[0] + r.d[0] * t2, r.o[1] + r.d[1] * t2, r.o[2] + r.d[2] * t2 };
        if (cone_ref_sdf(k, p1) < cone_ref_sdf(k, p2)) hi = t2; else lo = t1;
    }
    double t = 0.5 * (lo + hi), p[3] = { r.o[0] + r.d[0] * t, r.o[1] + r.d[1] * t, r.o[2] + r.d[2] * t };
    return cone_ref_sdf(k, p) / (px * fmax(-p[2], 1e-6));
}

/* How far p is from a flat cone's rims (where its normal jumps), world units */
static double cone_ref_rim_dist(const cone_ref_t* k, const double p[3]) {
    if (!k->flat) return 1e30;
    double ab[3] = { k->b[0] - k->a[0], k->b[1] - k->a[1], k->b[2] - k->a[2] }, L = sqrt(dot3(ab, ab));
    double u[3] = { ab[0] / L, ab[1] / L, ab[2] / L };
    double pa[3] = { p[0] - k->a[0], p[1] - k->a[1], p[2] - k->a[2] }, y = dot3(pa, u);
    double q[3] = { pa[0] - u[0] * y, pa[1] - u[1] * y, pa[2] - u[2] * y }, x = sqrt(dot3(q, q));
    return fmin(hypot(y, x - k->ra), hypot(y - L, x - k->rb));
}

UTEST(render, capsules_shaped_match_their_cones) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 192, H = 144;
    cam_t c = { .tan_half = 0.4f, .aspect = (float)W / H, .n = 0.1f };
    md_render_camera_t cam = make_camera(&c);
    const double pxd = 2.0 * c.tan_half / H;           /* a pixel per unit depth */

    /* Cylinders (a second instance turned and scaled 1.3), cones pointing
       across, away and at the eye, a truncated cone, round tapered cones,
       lines in pixels (round and flat ends; one long and receding, one
       through the near plane and behind the eye) */
    enum { NB = 6 };
    static const float pos[NB][3][2][3] = {
        { { { -5.0f, 3.0f, -14.0f }, { -1.0f, 3.2f, -14.0f } },         /* flat cylinders */
          { {  0.5f, 3.0f, -13.0f }, {  4.5f, 1.5f, -16.0f } },
          { { -4.8f, 0.5f, -12.0f }, { -4.5f, 0.8f, -16.5f } } },        /* nearly along the view: an end faces it */
        { { { -2.5f, -0.5f, -14.0f }, { -0.5f, 0.3f, -14.0f } },        /* cones */
          { {  0.5f, -1.5f, -12.0f }, {  0.8f, -1.2f, -15.0f } },        /* pointing away: its base faces the eye */
          { {  2.6f, -0.6f, -15.0f }, {  2.1f, -0.1f, -12.0f } } },      /* pointing at the eye */
        { { { -5.0f, -3.0f, -14.0f }, { -2.0f, -2.5f, -15.0f } },       /* a truncated cone */
          { { -1.5f, -3.2f, -13.0f }, { -1.2f, -3.0f, -16.0f } },
          { { -3.5f, 4.0f, -18.0f }, { -2.0f, 4.5f, -18.0f } } },
        { { {  0.0f, -3.0f, -13.0f }, {  3.0f, -3.5f, -15.0f } },       /* round, tapered */
          { {  3.5f,  2.5f, -12.0f }, {  3.8f,  2.8f, -16.0f } },
          { {  4.5f, -1.0f, -14.0f }, {  4.9f, -1.1f, -14.2f } } },      /* short: the cone between the balls steep */
        { { { -5.5f, -4.5f, -20.0f }, {  5.5f, -4.2f, -30.0f } },       /* in pixels, round ends */
          { {  1.0f, -0.5f, -20.0f }, {  0.3f, -0.2f,   2.0f } },        /* through the near plane */
          { { -6.0f,  5.0f, -60.0f }, { 60.0f, 9.0f, -80.0f } } },       /* long, thin, far */
        { { { -6.0f,  4.2f, -25.0f }, {  6.0f,  4.5f, -18.0f } },       /* in pixels, flat ends */
          { { -1.0f,  1.5f, -10.0f }, { -1.0f,  1.5f, -30.0f } },        /* along the view */
          { {  2.0f,  1.0f, -11.0f }, {  2.8f,  1.4f, -11.5f } } },
    };
    static const struct { bool flat, tapered; float radius, end_radius, px; } st[NB] = {
        { true,  false, 0.35f, 0.0f,  0.0f },
        { true,  true,  0.6f,  0.0f,  0.0f },
        { true,  true,  0.5f,  0.2f,  0.0f },
        { false, true,  0.45f, 0.15f, 0.0f },
        { false, false, 0.0f,  0.0f,  1.5f },
        { true,  false, 0.0f,  0.0f,  1.25f },
    };
    float xf[32] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };
    {   /* about (0, 0, -14): turned 30 degrees about y, scaled 1.3, moved back and down */
        const double ang = 30 * 3.14159265358979 / 180, sc = 1.3, ctr[3] = { 0, 0, -14 }, mv[3] = { 0.5, -1.0, -6.0 };
        float m[9] = { (float)(sc * cos(ang)), 0, (float)(-sc * sin(ang)), 0, (float)sc, 0, (float)(sc * sin(ang)), 0, (float)(sc * cos(ang)) };
        for (int i = 0; i < 9; ++i) xf[16 + 4 * (i / 3) + i % 3] = m[i];
        for (int i = 0; i < 3; ++i) xf[28 + i] = (float)(ctr[i] + mv[i] - (m[i] * ctr[0] + m[3 + i] * ctr[1] + m[6 + i] * ctr[2]));
    }
    md_render_batch_t b[NB];
    cone_ref_t ref[NB * 3 * 2];
    int nref = 0;
    for (int k = 0; k < NB; ++k) {
        memset(&b[k], 0, sizeof(b[k]));
        b[k].key = 60u + (uint32_t)k; b[k].kind = MD_RENDER_CAPSULES; b[k].count = 3; b[k].points.count = 6;
        b[k].points.position = f_buf(&fx, pos[k], sizeof(pos[k]));
        b[k].style.capsules.radius     = st[k].radius;
        b[k].style.capsules.flat       = st[k].flat;
        b[k].style.capsules.tapered    = st[k].tapered;
        b[k].style.capsules.end_radius = st[k].end_radius;
        b[k].style.capsules.radius_px  = st[k].px;
        if (k == 0) { b[k].transform = xf; b[k].transform_count = 2; }
        for (int in = 0; in < (k == 0 ? 2 : 1); ++in) {
            const float* m = xf + 16 * in;
            const double s = sqrt((double)m[0] * m[0] + (double)m[1] * m[1] + (double)m[2] * m[2]);
            for (int e = 0; e < 3; ++e) {
                cone_ref_t* q = &ref[nref++];
                for (int j = 0; j < 2; ++j) for (int i = 0; i < 3; ++i) {
                    const float* v = pos[k][e][j];
                    (j ? q->b : q->a)[i] = (double)m[i] * v[0] + (double)m[4 + i] * v[1] + (double)m[8 + i] * v[2] + m[12 + i];
                }
                q->flat = st[k].flat;
                q->key = b[k].key; q->linear = (uint32_t)(in * 3 + e);
                if (st[k].px > 0) {
                    /* cut at the near plane, cap_px pixels at each end's depth */
                    double zn = -c.n;
                    if (q->a[2] > zn) { double t = (zn - q->b[2]) / (q->a[2] - q->b[2]); for (int i = 0; i < 3; ++i) q->a[i] = q->b[i] + (q->a[i] - q->b[i]) * t; }
                    if (q->b[2] > zn) { double t = (zn - q->a[2]) / (q->b[2] - q->a[2]); for (int i = 0; i < 3; ++i) q->b[i] = q->a[i] + (q->b[i] - q->a[i]) * t; }
                    q->ra = st[k].px * pxd * -q->a[2] * s;
                    q->rb = st[k].px * pxd * -q->b[2] * s;
                } else {
                    q->ra = st[k].radius * s;
                    q->rb = (st[k].tapered ? st[k].end_radius : st[k].radius) * s;
                }
            }
        }
    }
    md_render_frame_t f = frame_of(W, H, cam, b, NB);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));

    md_render_targets_t tg = md_render_view_targets(fx.view);
    uint32_t* id = malloc(W * H * 8);
    float*    dp = malloc(W * H * 4);
    uint32_t* nb = malloc(W * H * 4);
    ASSERT_TRUE(f_read(&fx, tg.id, id, W * H * 8) && f_read(&fx, tg.depth, dp, W * H * 4) && f_read(&fx, tg.normal, nb, W * H * 4));
    int bad = 0, covered = 0, judged = 0, nerr = 0, ends = 0;
    int per_key[NB] = {0};
    double worst = 0, margin_worst = 0;
    double* errs = malloc(W * H * sizeof(double));
    const double M = 0.1;                  /* pixels: the coverage judged this far from the silhouette */
    for (uint32_t y = 0; y < H; ++y) {
        for (uint32_t x = 0; x < W; ++x) {
            uint32_t i = y * W + x;
            ray_t r = pixel_ray(&c, W, H, x + 0.5, y + 0.5);
            double mmin = 1e30;
            for (int k = 0; k < nref; ++k) mmin = fmin(mmin, cone_ref_ray_min(&ref[k], r, pxd));
            uint32_t iy = id[2 * i + 1];
            if (iy == 0) {
                if (mmin < -M) { bad++; if (bad < 8) printf("    pixel (%u, %u): empty, %.3f px inside a cone\n", x, y, -mmin); }
                continue;
            }
            covered++;
            if (mmin > M) { bad++; if (bad < 8) printf("    pixel (%u, %u): covered, %.3f px off every cone\n", x, y, mmin); continue; }
            if (mmin > 0) margin_worst = fmax(margin_worst, mmin);
            double zv = -c.n / dp[i], tp = zv / r.d[2];
            double P[3] = { r.o[0] + r.d[0] * tp, r.o[1] + r.d[1] * tp, r.o[2] + r.d[2] * tp };
            double pxw = -P[2] * pxd;
            int who = -1;
            double best = 1e30, second = 1e30;
            for (int k = 0; k < nref; ++k) {
                double sd = cone_ref_sdf(&ref[k], P);
                if (sd < best) { second = best; best = sd; who = k; } else if (sd < second) second = sd;
            }
            worst = fmax(worst, fabs(best) / pxw);
            if (fabs(best) > 0.05 * pxw) { bad++; if (bad < 8) printf("    pixel (%u, %u): %.3f px off the cones\n", x, y, best / pxw); continue; }
            if (second - best > pxw) {
                judged++;
                uint32_t key = (iy >> 20) - 1u;     /* the slot: batches in order */
                if (id[2 * i] != ref[who].linear || key != ref[who].key - 60u) {
                    bad++; if (bad < 8) printf("    pixel (%u, %u): batch %u element %u, expected %u, %u\n", x, y, key, id[2 * i], ref[who].key - 60u, ref[who].linear);
                    continue;
                }
                per_key[key]++;
                /* the normal: the distance's gradient, away from a flat end's rim */
                double g[3], ng[3];
                for (int k = 0; k < 3; ++k) {
                    double a[3] = { P[0], P[1], P[2] }, bq[3] = { P[0], P[1], P[2] };
                    a[k] += 1e-6; bq[k] -= 1e-6;
                    g[k] = (cone_ref_sdf(&ref[who], a) - cone_ref_sdf(&ref[who], bq)) / 2e-6;
                }
                double gl = sqrt(dot3(g, g));
                if (cone_ref_rim_dist(&ref[who], P) > 0.7 * pxw && fabs(dot3(g, r.d)) > 0.3 * gl) {
                    oct_decode(nb[i], ng);
                    double cs = dot3(g, ng) / gl;
                    errs[nerr++] = acos(cs > 1 ? 1 : cs < -1 ? -1 : cs) * 180.0 / 3.14159265358979;
                    double ab[3] = { ref[who].b[0] - ref[who].a[0], ref[who].b[1] - ref[who].a[1], ref[who].b[2] - ref[who].a[2] };
                    ends += ref[who].flat && fabs(dot3(g, ab)) > 0.999 * gl * sqrt(dot3(ab, ab));
                }
            }
        }
    }
    qsort(errs, (size_t)nerr, sizeof(double), cmp_double);
    double p50 = nerr ? errs[nerr / 2] : 0, p99 = nerr ? errs[(nerr * 99) / 100] : 0;
    printf("    %d pixels covered, %d clear (by batch %d %d %d %d %d %d); the farthest %.4f px from the cones, %.4f px outside; normals over %d (%d on flat ends): median %.3f, 99th %.3f degrees\n",
           covered, judged, per_key[0], per_key[1], per_key[2], per_key[3], per_key[4], per_key[5], worst, margin_worst, nerr, ends, p50, p99);
    EXPECT_EQ(0, bad);
    EXPECT_GT(covered, 2000);
    EXPECT_GT(judged, covered * 3 / 4);
    for (int k = 0; k < NB; ++k) EXPECT_GT(per_key[k], 40);
    EXPECT_GT(ends, 30);
    EXPECT_LT(p50, 0.2);
    EXPECT_LT(p99, 1.0);
    free(errs);

    /* Picks: near each end of the first cylinder (its nearer point's
       attribute), and on the long line in pixels */
    for (int end = 0; end < 2; ++end) {
        const cone_ref_t* q = &ref[0];
        double s = end ? 0.85 : 0.15, p[3];
        for (int k = 0; k < 3; ++k) p[k] = q->a[k] + (q->b[k] - q->a[k]) * s;
        uint32_t px = (uint32_t)((p[0] / -p[2] / (c.tan_half * c.aspect) * 0.5 + 0.5) * W);
        uint32_t py = (uint32_t)((0.5 - p[1] / -p[2] / c.tan_half * 0.5) * H);
        md_render_pick_t pk;
        ASSERT_TRUE(f_pick(&fx, px, py, &pk));
        EXPECT_EQ(60u, (unsigned)pk.key);
        EXPECT_EQ(0u, pk.element);
        EXPECT_EQ((uint32_t)end, pk.attr);
    }
    {
        const cone_ref_t* q = &ref[3 * 2 + 3 * 3];      /* batch 4's first */
        double p[3];
        for (int k = 0; k < 3; ++k) p[k] = q->a[k] + (q->b[k] - q->a[k]) * 0.7;
        double sx = (p[0] / -p[2] / (c.tan_half * c.aspect) * 0.5 + 0.5) * W, sy = (0.5 - p[1] / -p[2] / c.tan_half * 0.5) * H;
        md_render_pick_t pk;
        ASSERT_TRUE(f_pick(&fx, (uint32_t)sx, (uint32_t)sy, &pk));
        EXPECT_EQ(64u, (unsigned)pk.key);
        EXPECT_EQ(0u, pk.element);
    }

    /* Refused: shaped and dashed, tapered in pixels, a negative radius in pixels */
    for (int k = 0; k < 3; ++k) {
        md_render_batch_t x = b[k == 0 ? 0 : k == 1 ? 4 : 5];
        if (k == 0) { x.style.capsules.dash_length = 0.3f; x.style.capsules.gap_length = 0.2f; }
        if (k == 1) { x.style.capsules.tapered = true; x.style.capsules.end_radius = 0.1f; }
        if (k == 2) x.style.capsules.radius_px = -1.0f;
        md_render_frame_t fr = frame_of(W, H, cam, &x, 1);
        EXPECT_FALSE(md_render_view_draw(fx.gfx, fx.view, &fr, NULL));
    }
    free(id); free(dp); free(nb);
    f_close(&fx);
}

/* ---- Triangles ------------------------------------------------------------------ */

/* An icosphere: the icosahedron's faces split in four `levels` times, the new
   corners pushed out onto the unit sphere; counter-clockwise from outside */
typedef struct { float* xyz; uint32_t nv, cap_v; uint32_t* tri; uint32_t nt; } ico_t;

static uint32_t ico_mid(ico_t* m, uint32_t a, uint32_t b) {
    float p[3];
    for (int k = 0; k < 3; ++k) p[k] = 0.5f * (m->xyz[3 * a + k] + m->xyz[3 * b + k]);
    float l = sqrtf(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    for (int k = 0; k < 3; ++k) p[k] /= l;
    for (uint32_t i = 0; i < m->nv; ++i) {
        if (fabsf(m->xyz[3 * i] - p[0]) < 1e-6f && fabsf(m->xyz[3 * i + 1] - p[1]) < 1e-6f && fabsf(m->xyz[3 * i + 2] - p[2]) < 1e-6f) return i;
    }
    for (int k = 0; k < 3; ++k) m->xyz[3 * m->nv + k] = p[k];
    return m->nv++;
}

static ico_t ico_make(int levels) {
    const float t = (1.0f + sqrtf(5.0f)) / 2.0f;
    static const float base[12][3] = { {-1, 0, 0}, {1, 0, 0}, {-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1}, {0, 0, -1}, {0, 0, 1} };
    static const uint32_t faces[20][3] = { {0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11}, {1, 5, 9}, {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
                                           {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8}, {3, 8, 9}, {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1} };
    (void)base;
    ico_t m;
    uint32_t nt = 20;
    for (int l = 0; l < levels; ++l) nt *= 4;
    m.cap_v = nt * 3;
    m.xyz = malloc(m.cap_v * 3 * sizeof(float));
    m.tri = malloc(nt * 3 * sizeof(uint32_t));
    const float v12[12][3] = { {-1, t, 0}, {1, t, 0}, {-1, -t, 0}, {1, -t, 0}, {0, -1, t}, {0, 1, t}, {0, -1, -t}, {0, 1, -t}, {t, 0, -1}, {t, 0, 1}, {-t, 0, -1}, {-t, 0, 1} };
    m.nv = 12;
    for (int i = 0; i < 12; ++i) {
        float l = sqrtf(v12[i][0] * v12[i][0] + v12[i][1] * v12[i][1] + v12[i][2] * v12[i][2]);
        for (int k = 0; k < 3; ++k) m.xyz[3 * i + k] = v12[i][k] / l;
    }
    memcpy(m.tri, faces, sizeof(faces));
    m.nt = 20;
    for (int l = 0; l < levels; ++l) {
        uint32_t* nt2 = malloc(m.nt * 4 * 3 * sizeof(uint32_t));
        for (uint32_t f = 0; f < m.nt; ++f) {
            uint32_t a = m.tri[3 * f], b = m.tri[3 * f + 1], c = m.tri[3 * f + 2];
            uint32_t ab = ico_mid(&m, a, b), bc = ico_mid(&m, b, c), ca = ico_mid(&m, c, a);
            uint32_t q[12] = { a, ab, ca, b, bc, ab, c, ca, bc, ab, bc, ca };
            memcpy(nt2 + 12 * f, q, sizeof(q));
        }
        free(m.tri);
        m.tri = nt2;
        m.nt *= 4;
    }
    return m;
}

typedef struct {
    double   q[3][3];                  /* view space (the camera is the identity) */
    double   n[3][3];                  /* the corners' normals; flat: none */
    bool     smooth;
    uint32_t col[3];                   /* the corners' colours */
    float    sharp;
    uint32_t slot, linear, attr[3];
} tri_ref_t;

/* Moller-Trumbore, double: t > 0 and the weights (1 - u - v, u, v) */
static bool tri_ref_hit(const tri_ref_t* k, ray_t r, double* t, double w[3]) {
    double e1[3], e2[3], pv[3], tv[3], qv[3];
    for (int i = 0; i < 3; ++i) { e1[i] = k->q[1][i] - k->q[0][i]; e2[i] = k->q[2][i] - k->q[0][i]; tv[i] = r.o[i] - k->q[0][i]; }
    pv[0] = r.d[1] * e2[2] - r.d[2] * e2[1]; pv[1] = r.d[2] * e2[0] - r.d[0] * e2[2]; pv[2] = r.d[0] * e2[1] - r.d[1] * e2[0];
    double det = dot3(e1, pv);
    if (fabs(det) < 1e-14) return false;
    double u = dot3(tv, pv) / det;
    qv[0] = tv[1] * e1[2] - tv[2] * e1[1]; qv[1] = tv[2] * e1[0] - tv[0] * e1[2]; qv[2] = tv[0] * e1[1] - tv[1] * e1[0];
    double v = dot3(r.d, qv) / det;
    if (u < 0 || v < 0 || u + v > 1) return false;
    *t = dot3(e2, qv) / det;
    w[0] = 1 - u - v; w[1] = u; w[2] = v;
    return *t > 0;
}

/* The nearest hit along a ray: its triangle (-1 none), t and weights */
static int tri_ref_cast(const tri_ref_t* refs, int n, ray_t r, double* t, double w[3]) {
    int best = -1;
    double tb = 1e30, wb[3] = {0};
    for (int k = 0; k < n; ++k) {
        double tk, wk[3];
        if (tri_ref_hit(&refs[k], r, &tk, wk) && tk < tb) { tb = tk; best = k; memcpy(wb, wk, sizeof(wb)); }
    }
    *t = tb;
    memcpy(w, wb, sizeof(wb));
    return best;
}

UTEST(render, triangles_match_their_mesh) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 192, H = 144;
    cam_t c = { .tan_half = 0.4f, .aspect = (float)W / H, .n = 0.1f };
    md_render_camera_t cam = make_camera(&c);

    /* Batch 0: an icosphere, its points and normals float4 (as md_gto_iso's
       meshes), indexed, smooth, marked (state 0 on every corner); batch 1:
       the back half of another, open toward the eye (its inside, back faces
       with outward normals, turned to the eye); batch 2: a flat tetrahedron
       of corners coloured red, green and blue, sharp (the nearest corner's),
       twice (the second instance turned and moved); batch 3: the same
       blended (sharpness 0); batch 4: a floor of two triangles reaching
       behind the eye (clipped at the near plane) */
    ico_t s = ico_make(2);
    const float c0[3] = { -2.6f, 1.0f, -14.0f }, r0 = 2.0f;
    const float c1[3] = {  2.8f, 1.2f, -13.0f }, r1 = 1.8f;
    float* pos0 = malloc(s.nv * 16); float* nrm0 = malloc(s.nv * 16);
    float* pos1 = malloc(s.nv * 12); float* nrm1 = malloc(s.nv * 12);
    for (uint32_t i = 0; i < s.nv; ++i) {
        for (int k = 0; k < 3; ++k) {
            pos0[4 * i + k] = c0[k] + r0 * s.xyz[3 * i + k]; nrm0[4 * i + k] = s.xyz[3 * i + k];
            pos1[3 * i + k] = c1[k] + r1 * s.xyz[3 * i + k]; nrm1[3 * i + k] = s.xyz[3 * i + k];
        }
        pos0[4 * i + 3] = 1.0f; nrm0[4 * i + 3] = 0.25f;     /* w unread */
    }
    uint32_t* back = malloc(s.nt * 12);
    uint32_t nback = 0;
    for (uint32_t f = 0; f < s.nt; ++f) {
        bool behind = true;
        for (int k = 0; k < 3; ++k) behind = behind && s.xyz[3 * s.tri[3 * f + k] + 2] <= 0.05f;
        if (behind) { memcpy(back + 3 * nback, s.tri + 3 * f, 12); ++nback; }
    }
    const float tet[4][3] = { { -1.6f, -2.6f, -11.0f }, { 0.4f, -2.2f, -12.5f }, { -0.8f, -0.6f, -12.0f }, { -0.4f, -2.0f, -10.0f } };
    const uint32_t tet_idx[12] = { 0, 1, 2,  0, 3, 1,  1, 3, 2,  2, 3, 0 };
    const uint32_t tet_col[4] = { 0xFF2020E0u, 0xFF20E020u, 0xFFE02020u, 0xFF20E0E0u };
    const float floor_[4][3] = { { -6.0f, -4.0f, 3.0f }, { 6.0f, -4.0f, 3.0f }, { 6.0f, -4.0f, -40.0f }, { -6.0f, -4.0f, -40.0f } };
    const uint32_t floor_idx[6] = { 0, 1, 2,  0, 2, 3 };
    float xf[32] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };
    {   /* turned 40 degrees about z through the tetrahedron, moved right and back */
        const double a = 40 * 3.14159265358979 / 180, ctr[3] = { -0.6, -1.85, -11.4 }, mv[3] = { 3.2, -0.2, -1.5 };
        float m[9] = { (float)cos(a), (float)sin(a), 0, (float)-sin(a), (float)cos(a), 0, 0, 0, 1 };
        for (int i = 0; i < 9; ++i) xf[16 + 4 * (i / 3) + i % 3] = m[i];
        for (int i = 0; i < 3; ++i) xf[28 + i] = (float)(ctr[i] + mv[i] - (m[i] * ctr[0] + m[3 + i] * ctr[1] + m[6 + i] * ctr[2]));
    }
    uint8_t* marks0 = malloc(s.nv);
    memset(marks0, MD_RENDER_MARK_STATE0, s.nv);

    enum { NB = 5 };
    md_render_batch_t b[NB];
    memset(b, 0, sizeof(b));
    for (int k = 0; k < NB; ++k) { b[k].key = 80u + (uint32_t)k; b[k].kind = MD_RENDER_TRIANGLES; b[k].color.mode = MD_RENDER_COLOR_UNIFORM; }
    b[0].count = s.nt; b[0].points.count = s.nv;
    b[0].points.position = f_buf(&fx, pos0, s.nv * 16); b[0].points.stride = 16;
    b[0].index = f_buf(&fx, s.tri, s.nt * 12);
    b[0].style.triangles.normal = f_buf(&fx, nrm0, s.nv * 16); b[0].style.triangles.normal_stride = 16;
    b[0].points.marks = f_buf(&fx, marks0, s.nv);
    b[0].color.uniform = 0xFF3070D0u;
    b[1].count = nback; b[1].points.count = s.nv;
    b[1].points.position = f_buf(&fx, pos1, s.nv * 12);
    b[1].index = f_buf(&fx, back, nback * 12);
    b[1].style.triangles.normal = f_buf(&fx, nrm1, s.nv * 12);
    b[1].color.uniform = 0xFF40C0F0u;
    md_gpu_addr_t TP = f_buf(&fx, tet, sizeof(tet)), TI = f_buf(&fx, tet_idx, sizeof(tet_idx)), TC = f_buf(&fx, tet_col, sizeof(tet_col));
    for (int k = 2; k < 4; ++k) {
        b[k].count = 4; b[k].points.count = 4; b[k].points.position = TP; b[k].index = TI;
        b[k].points.color = TC; b[k].color.mode = MD_RENDER_COLOR_POINTS;
    }
    b[2].color.sharpness = 1.0f; b[2].transform = xf; b[2].transform_count = 2;
    b[3].color.sharpness = 0.0f;
    {   /* the blended one beside the sharp one's first instance */
        static float xf3[16] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,3.4f,0,1 };
        b[3].transform = xf3; b[3].transform_count = 1;
    }
    b[4].count = 2; b[4].points.count = 4; b[4].points.position = f_buf(&fx, floor_, sizeof(floor_)); b[4].index = f_buf(&fx, floor_idx, sizeof(floor_idx));
    b[4].color.uniform = 0xFF909090u;

    /* The references, view space (the camera at the origin looking down -z) */
    int nref = 0, cap = (int)(s.nt + nback + 4 * 3 + 2);
    tri_ref_t* ref = calloc((size_t)cap, sizeof(tri_ref_t));
    for (int k = 0; k < NB; ++k) {
        const uint32_t ninst = b[k].transform_count ? b[k].transform_count : 1;
        for (uint32_t in = 0; in < ninst; ++in) {
            const float* m = b[k].transform ? b[k].transform + 16 * in : xf;
            for (uint32_t e = 0; e < b[k].count; ++e) {
                tri_ref_t* q = &ref[nref++];
                q->slot = (uint32_t)k; q->linear = in * b[k].count + e;
                q->smooth = k <= 1; q->sharp = b[k].color.sharpness;
                for (int j = 0; j < 3; ++j) {
                    uint32_t pi = k == 0 ? s.tri[3 * e + j] : k == 1 ? back[3 * e + j] : k <= 3 ? tet_idx[3 * e + j] : floor_idx[3 * e + j];
                    const float* pp = k == 0 ? pos0 + 4 * pi : k == 1 ? pos1 + 3 * pi : k <= 3 ? tet[pi] : floor_[pi];
                    for (int i = 0; i < 3; ++i) q->q[j][i] = (double)m[i] * pp[0] + (double)m[4 + i] * pp[1] + (double)m[8 + i] * pp[2] + m[12 + i];
                    if (k <= 1) for (int i = 0; i < 3; ++i) q->n[j][i] = s.xyz[3 * pi + i];
                    q->col[j]  = k == 0 ? b[0].color.uniform : k == 1 ? b[1].color.uniform : k <= 3 ? tet_col[pi] : b[4].color.uniform;
                    q->attr[j] = pi;
                }
            }
        }
    }

    md_render_frame_t f = frame_of(W, H, cam, b, NB);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    md_render_targets_t tg = md_render_view_targets(fx.view);
    uint32_t* id = malloc(W * H * 8);
    float*    dp = malloc(W * H * 4);
    uint32_t* nb = malloc(W * H * 4);
    uint32_t* cl = malloc(W * H * 4);
    uint32_t* mk = malloc(W * H * 4);
    ASSERT_TRUE(f_read(&fx, tg.id, id, W * H * 8) && f_read(&fx, tg.depth, dp, W * H * 4) && f_read(&fx, tg.normal, nb, W * H * 4) &&
                f_read(&fx, tg.color, cl, W * H * 4) && f_read(&fx, tg.mask, mk, W * H * 4));

    int bad = 0, sure = 0, per_slot[NB] = {0}, nerr = 0, backs = 0, marked = 0, sharp_ok = 0, blend_ok = 0;
    double worst_depth = 0, worst_col = 0;
    double* errs = malloc(W * H * sizeof(double));
    const double pxd = 2.0 * c.tan_half / H;
    for (uint32_t y = 0; y < H; ++y) {
        for (uint32_t x = 0; x < W; ++x) {
            const uint32_t i = y * W + x;
            /* sure: the centre's triangle (or none) is that of four rays 0.05 px about it too */
            double t, w[3];
            ray_t r = pixel_ray(&c, W, H, x + 0.5, y + 0.5);
            int k0 = tri_ref_cast(ref, nref, r, &t, w);
            bool agree = true;
            static const double off[4][2] = { { 0.05, 0 }, { -0.05, 0 }, { 0, 0.05 }, { 0, -0.05 } };
            for (int o = 0; o < 4 && agree; ++o) {
                double to, wo[3];
                agree = tri_ref_cast(ref, nref, pixel_ray(&c, W, H, x + 0.5 + off[o][0], y + 0.5 + off[o][1]), &to, wo) == k0;
            }
            if (!agree) continue;
            sure++;
            const uint32_t iy = id[2 * i + 1];
            if (k0 < 0) {
                if (iy != 0) { bad++; if (bad < 8) printf("    pixel (%u, %u): covered, no triangle there\n", x, y); }
                continue;
            }
            const tri_ref_t* q = &ref[k0];
            if (iy == 0 || (iy >> 20) - 1u != q->slot || id[2 * i] != q->linear) {
                bad++; if (bad < 8) printf("    pixel (%u, %u): batch %d element %u, expected %u, %u\n", x, y, iy ? (int)(iy >> 20) - 1 : -1, id[2 * i], q->slot, q->linear);
                continue;
            }
            per_slot[q->slot]++;
            /* depth: the raster's own, within a hundredth of a pixel's depth step at it */
            const double zv = -c.n / dp[i], tp = zv / r.d[2];
            const double dz = fabs(tp - t) * fabs(r.d[2]);
            const double tol = 1e-5 * -zv;                             /* far more than a float depth's step there */
            worst_depth = fmax(worst_depth, dz / tol);
            if (dz > tol) { bad++; if (bad < 8) printf("    pixel (%u, %u): depth %.6f against %.6f\n", x, y, -zv, t * -r.d[2]); continue; }
            /* the normal: by the weights (smooth) or the plane's, facing the eye */
            double n[3] = {0, 0, 0};
            if (q->smooth) {
                for (int j = 0; j < 3; ++j) for (int a = 0; a < 3; ++a) n[a] += w[j] * q->n[j][a];
            } else {
                double e1[3], e2[3];
                for (int a = 0; a < 3; ++a) { e1[a] = q->q[1][a] - q->q[0][a]; e2[a] = q->q[2][a] - q->q[0][a]; }
                n[0] = e1[1] * e2[2] - e1[2] * e2[1]; n[1] = e1[2] * e2[0] - e1[0] * e2[2]; n[2] = e1[0] * e2[1] - e1[1] * e2[0];
            }
            double nl = sqrt(dot3(n, n));
            for (int a = 0; a < 3; ++a) n[a] /= nl;
            if (dot3(n, r.d) > 0) { for (int a = 0; a < 3; ++a) n[a] = -n[a]; backs += q->slot == 1; }
            double ng[3];
            oct_decode(nb[i], ng);
            double cs = dot3(n, ng);
            errs[nerr++] = acos(cs > 1 ? 1 : cs < -1 ? -1 : cs) * 180.0 / 3.14159265358979;
            /* the colour: the nearest corner's (sharp, clear of a tie) or the blend */
            double ws[3] = { w[0], w[1], w[2] };
            int km = w[1] > w[0] ? (w[2] > w[1] ? 2 : 1) : (w[2] > w[0] ? 2 : 0);
            double second = 0;
            for (int j = 0; j < 3; ++j) if (j != km) second = fmax(second, w[j]);
            if (q->sharp >= 1.0f) {
                if (w[km] - second < 0.03) continue;
                for (int j = 0; j < 3; ++j) ws[j] = j == km ? 1 : 0;
            }
            double ce = 0;
            for (int ch = 0; ch < 3; ++ch) {
                double e = 0;
                for (int j = 0; j < 3; ++j) e += ws[j] * ((q->col[j] >> (8 * ch)) & 0xFF);
                ce = fmax(ce, fabs(e - ((cl[i] >> (8 * ch)) & 0xFF)));
            }
            worst_col = fmax(worst_col, ce);
            if (ce > 2.5) { bad++; if (bad < 8) printf("    pixel (%u, %u): colour %08x, expected another (%.1f off)\n", x, y, cl[i], ce); continue; }
            sharp_ok += q->sharp >= 1.0f;
            blend_ok += q->sharp < 1.0f && q->slot == 3;
            /* the mask: batch 0's every corner in state 0, visible */
            if (q->slot == 0) marked += (mk[i] & 0xFF) >= 128;
        }
    }
    qsort(errs, (size_t)nerr, sizeof(double), cmp_double);
    const double p50 = nerr ? errs[nerr / 2] : 0, p99 = nerr ? errs[(nerr * 99) / 100] : 0;
    printf("    %d pixels sure (by batch %d %d %d %d %d), %d of the bowl's from behind; depth at most %.3f of its tolerance, colour %.2f off; normals median %.4f, 99th %.4f degrees; %d sharp, %d blended, %d marked\n",
           sure, per_slot[0], per_slot[1], per_slot[2], per_slot[3], per_slot[4], backs, worst_depth, worst_col, p50, p99, sharp_ok, blend_ok, marked);
    EXPECT_EQ(0, bad);
    for (int k = 0; k < NB; ++k) EXPECT_GT(per_slot[k], 100);
    EXPECT_GT(backs, 100);
    EXPECT_GT(sharp_ok, 100);
    EXPECT_GT(blend_ok, 100);
    EXPECT_EQ(per_slot[0], marked);
    EXPECT_LT(p50, 0.05);
    EXPECT_LT(p99, 0.5);
    free(errs);

    /* A pick on the tetrahedron's second instance: its batch, instance, face
       and the corner of the largest weight */
    bool picked = false;
    for (int k = 0; k < nref; ++k) {
        if (ref[k].slot != 2 || ref[k].linear < 4) continue;
        double ctr[3];
        for (int a = 0; a < 3; ++a) ctr[a] = (0.6 * ref[k].q[0][a] + 0.2 * ref[k].q[1][a] + 0.2 * ref[k].q[2][a]);
        const uint32_t px = (uint32_t)((ctr[0] / -ctr[2] / (c.tan_half * c.aspect) * 0.5 + 0.5) * W);
        const uint32_t py = (uint32_t)((0.5 - ctr[1] / -ctr[2] / c.tan_half * 0.5) * H);
        double t, w[3];
        int kh = tri_ref_cast(ref, nref, pixel_ray(&c, W, H, px + 0.5, py + 0.5), &t, w);
        if (kh != k) continue;      /* that face hidden there: the next */
        md_render_pick_t pk;
        ASSERT_TRUE(f_pick(&fx, px, py, &pk));
        int km = w[1] > w[0] ? (w[2] > w[1] ? 2 : 1) : (w[2] > w[0] ? 2 : 0);
        EXPECT_EQ(82u, (unsigned)pk.key);
        EXPECT_EQ(1u, pk.instance);
        EXPECT_EQ(ref[k].linear - 4u, pk.element);
        EXPECT_EQ(ref[k].attr[km], pk.attr);
        EXPECT_NEAR(t * -pixel_ray(&c, W, H, px + 0.5, py + 0.5).d[2], (double)-pk.position[2], 1e-3);
        picked = true;
        break;
    }
    EXPECT_TRUE(picked);

    /* Motion: the blended tetrahedron's corners given velocities of a
       turning about its centre and a drift, so the motion at a point is the
       field's there, through its corners by the weights. A raster at time 1
       after one at 0 (its positions as they were: the velocity alone says
       where it was) */
    {
        const double om[3] = { 0.3, -0.2, 0.5 }, v0[3] = { 0.4, 0.1, -0.3 }, cc[3] = { -0.6, -1.85, -11.4 };
        uint32_t pv[8];
        for (int j = 0; j < 4; ++j) {
            double d[3] = { tet[j][0] - cc[0], tet[j][1] - cc[1], tet[j][2] - cc[2] };
            float vj[3] = { (float)(om[1] * d[2] - om[2] * d[1] + v0[0]), (float)(om[2] * d[0] - om[0] * d[2] + v0[1]), (float)(om[0] * d[1] - om[1] * d[0] + v0[2]) };
            md_render_pack_velocity(vj, pv + 2 * j);
        }
        b[3].points.velocity = f_buf(&fx, pv, sizeof(pv));
        f = frame_of(W, H, cam, b, NB);
        f.time = 0.0;
        ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
        f.time = 1.0;
        f.content[MD_RENDER_CONTENT_GEOMETRY]++;
        ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
        uint32_t* mo = malloc(W * H * 4);
        ASSERT_TRUE(f_read(&fx, tg.motion, mo, W * H * 4) && f_read(&fx, tg.id, id, W * H * 8));
        int moved = 0, still = 0, mbad = 0;
        double mworst = 0;
        for (uint32_t y = 0; y < H; ++y) {
            for (uint32_t x = 0; x < W; ++x) {
                const uint32_t i = y * W + x;
                if (id[2 * i + 1] == 0) continue;
                double t, w[3];
                ray_t r = pixel_ray(&c, W, H, x + 0.5, y + 0.5);
                int k = tri_ref_cast(ref, nref, r, &t, w);
                if (k < 0 || (id[2 * i + 1] >> 20) - 1u != ref[k].slot || id[2 * i] != ref[k].linear) continue;
                const double mx = half_bits_to_float(mo[i] & 0xFFFFu), my = half_bits_to_float(mo[i] >> 16);
                if (ref[k].slot != 3) { still += mx == 0.0 && my == 0.0; continue; }
                /* the hit point, the field there (the batch's transform: a move up by 3.4) */
                double P[3] = { r.o[0] + r.d[0] * t, r.o[1] + r.d[1] * t, r.o[2] + r.d[2] * t };
                double d[3] = { P[0] - cc[0], P[1] - 3.4 - cc[1], P[2] - cc[2] };
                double v[3] = { om[1] * d[2] - om[2] * d[1] + v0[0], om[2] * d[0] - om[0] * d[2] + v0[1], om[0] * d[1] - om[1] * d[0] + v0[2] };
                double Q[3] = { P[0] - v[0], P[1] - v[1], P[2] - v[2] };
                double ex = ((P[0] / -P[2] - Q[0] / -Q[2]) / (c.tan_half * c.aspect)) * 0.5 * W;
                double ey = -((P[1] / -P[2] - Q[1] / -Q[2]) / c.tan_half) * 0.5 * H;
                double e = fmax(fabs(ex - mx), fabs(ey - my));
                mworst = fmax(mworst, e);
                if (e > 0.02 + 0.002 * fmax(fabs(ex), fabs(ey))) { mbad++; if (mbad < 6) printf("    pixel (%u, %u): motion (%.3f, %.3f), expected (%.3f, %.3f)\n", x, y, mx, my, ex, ey); }
                moved++;
            }
        }
        printf("    motion: %d pixels moving, at most %.4f px off; %d still\n", moved, mworst, still);
        EXPECT_EQ(0, mbad);
        EXPECT_GT(moved, 100);
        EXPECT_GT(still, 1000);
        free(mo);
    }

    /* Refused: a fade, a normal stride not of whole floats */
    for (int k = 0; k < 2; ++k) {
        md_render_batch_t x = b[0];
        if (k == 0) { x.fade.mode = MD_RENDER_FADE_SHRINK; x.fade.out[0] = 10.0f; x.fade.out[1] = 5.0f; }
        if (k == 1) x.style.triangles.normal_stride = 14;
        md_render_frame_t fr = frame_of(W, H, cam, &x, 1);
        EXPECT_FALSE(md_render_view_draw(fx.gfx, fx.view, &fr, NULL));
    }
    (void)pxd;
    free(id); free(dp); free(nb); free(cl); free(mk); free(ref);
    free(pos0); free(nrm0); free(pos1); free(nrm1); free(back); free(marks0); free(s.xyz); free(s.tri);
    f_close(&fx);
}

/* ---- Transparency ----------------------------------------------------------------- */

static bool f_read_buf(fix_t* f, md_gpu_addr_t src, void* dst, size_t size) {
    md_gpu_mem_t rb = md_gpu_malloc(f->gfx, MD_GPU_MEM_HOST_READ, size);
    if (!rb.cpu) return false;
    bool ok = md_gpu_copy(f->gfx, rb.gpu, src, size);
    md_gpu_stream_sync(f->gfx);
    if (ok) memcpy(dst, rb.cpu, size);
    md_gpu_free(f->gfx, rb.gpu);
    return ok;
}

/* A transparent primitive of the reference: a sphere, a capsule (round caps)
   or a triangle of smooth normals; its group and its batch's colour and
   opacity */
typedef struct {
    int      kind;                     /* 0 sphere, 1 capsule, 2 triangle */
    double   a[3], b[3], r;            /* sphere: a, r; capsule: a, b, r */
    tri_ref_t tri;
    uint32_t group, color;
    double   opacity;
} tp_ref_t;

/* Its nearest hit ahead, t and the normal there (outward); false: none */
static bool tp_hit(const tp_ref_t* q, ray_t r, double* t, double n[3]) {
    if (q->kind == 0) {
        double oc[3] = { r.o[0] - q->a[0], r.o[1] - q->a[1], r.o[2] - q->a[2] };
        double b = dot3(oc, r.d), c = dot3(oc, oc) - q->r * q->r, h = b * b - c;
        if (h < 0) return false;
        *t = -b - sqrt(h);
        if (*t <= 0) return false;
        for (int i = 0; i < 3; ++i) n[i] = (r.o[i] + r.d[i] * *t - q->a[i]) / q->r;
        return true;
    }
    if (q->kind == 1) {                /* Quilez's capsule */
        double ba[3] = { q->b[0] - q->a[0], q->b[1] - q->a[1], q->b[2] - q->a[2] };
        double oa[3] = { r.o[0] - q->a[0], r.o[1] - q->a[1], r.o[2] - q->a[2] };
        double baba = dot3(ba, ba), bard = dot3(ba, r.d), baoa = dot3(ba, oa), rdoa = dot3(r.d, oa), oaoa = dot3(oa, oa);
        double A = baba - bard * bard, B = baba * rdoa - baoa * bard, C = baba * oaoa - baoa * baoa - q->r * q->r * baba;
        double h = B * B - A * C;
        double tt = -1;
        if (h >= 0) {
            tt = (-B - sqrt(h)) / A;
            double y = baoa + tt * bard;
            if (y <= 0 || y >= baba) {
                double oc[3];
                for (int i = 0; i < 3; ++i) oc[i] = y <= 0 ? oa[i] : r.o[i] - q->b[i];
                double bb = dot3(r.d, oc), cc = dot3(oc, oc) - q->r * q->r, hh = bb * bb - cc;
                tt = hh > 0 ? -bb - sqrt(hh) : -1;
            }
        }
        if (tt <= 0) return false;
        *t = tt;
        double p[3] = { r.o[0] + r.d[0] * tt, r.o[1] + r.d[1] * tt, r.o[2] + r.d[2] * tt };
        double pa[3] = { p[0] - q->a[0], p[1] - q->a[1], p[2] - q->a[2] };
        double k = dot3(pa, ba) / baba;
        k = k < 0 ? 0 : k > 1 ? 1 : k;
        for (int i = 0; i < 3; ++i) n[i] = (pa[i] - ba[i] * k) / q->r;
        return true;
    }
    double w[3];
    if (!tri_ref_hit(&q->tri, r, t, w)) return false;
    for (int i = 0; i < 3; ++i) n[i] = w[0] * q->tri.n[0][i] + w[1] * q->tri.n[1][i] + w[2] * q->tri.n[2][i];
    double l = sqrt(dot3(n, n));
    for (int i = 0; i < 3; ++i) n[i] /= l;
    return true;
}

typedef struct { double t; int prim; } tp_frag_t;
static int tp_frag_cmp(const void* x, const void* y) {
    const tp_frag_t* a = x; const tp_frag_t* b = y;
    return a->t < b->t ? -1 : a->t > b->t ? 1 : 0;
}

/* The layers a ray should get: its fragments in front of the opaque one
   (t_op), the K nearest, then of each group the first. Out: prims, count. */
static int tp_layers(const tp_ref_t* q, int nq, ray_t r, double t_op, uint32_t K, int* out, double* out_t, int* nfrag) {
    tp_frag_t fr[512];
    int nf = 0;
    for (int i = 0; i < nq && nf < 512; ++i) {
        double t, n[3];
        if (tp_hit(&q[i], r, &t, n) && t < t_op) { fr[nf].t = t; fr[nf].prim = i; ++nf; }
    }
    qsort(fr, (size_t)nf, sizeof(fr[0]), tp_frag_cmp);
    *nfrag = nf;
    int nl = 0;
    for (int k = 0; k < nf && k < (int)K; ++k) {
        bool seen = false;
        for (int j = 0; j < nl; ++j) seen = seen || q[out[j]].group == q[fr[k].prim].group;
        if (seen) continue;
        out[nl] = fr[k].prim;
        out_t[nl] = fr[k].t;
        ++nl;
    }
    /* the K-th and the next too close to tell apart: -1 */
    for (int k = 1; k < nf; ++k) if (fabs(fr[k].t - fr[k - 1].t) < 1e-4 * fr[k].t && q[fr[k].prim].group != q[fr[k - 1].prim].group) return -1;
    return nl;
}

UTEST(render, transparency_layers) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    md_gpu_device_info_t info;
    md_gpu_device_info(fx.dev, &info);
    if (!info.atomic_int64 || !info.fragment_stores) UTEST_SKIP("no 64-bit atomics returning the old value in fragment stages (the K-buffer's)");
    const uint32_t W = 160, H = 120, K = 4;
    cam_t c = { .tan_half = 0.4f, .aspect = (float)W / H, .n = 0.1f };
    md_render_camera_t cam = make_camera(&c);

    /* Opaque: a large sphere behind everything. Group 77, transparency 0.4:
       a ghost of two spheres and the capsule between them, two batches of
       different colours. An icosphere of smooth normals, transparency 0.5,
       without a group (its front and its back both layers). Six single
       spheres without a group, transparency 0.3, along one line of sight:
       more surfaces there than K. */
    ico_t s = ico_make(2);
    const float mc[3] = { 0.6f, -0.6f, -14.0f }, mr = 1.6f;
    float* mpos = malloc(s.nv * 12); float* mnrm = malloc(s.nv * 12);
    for (uint32_t i = 0; i < s.nv; ++i) for (int k = 0; k < 3; ++k) { mpos[3 * i + k] = mc[k] + mr * s.xyz[3 * i + k]; mnrm[3 * i + k] = s.xyz[3 * i + k]; }
    float opq[4] = { 0.0f, 0.0f, -22.0f, 0.0f };
    float ga[6] = { -2.0f, 0.5f, -12.0f,   0.0f, 1.5f, -12.5f };
    float st[18];
    const double dir[3] = { 0.24, -0.16, -1.0 };
    for (int j = 0; j < 6; ++j) for (int k = 0; k < 3; ++k) st[3 * j + k] = (float)(dir[k] * (9.0 + 1.1 * j));
    enum { NB = 10 };
    md_render_batch_t b[NB];
    memset(b, 0, sizeof(b));
    for (int k = 0; k < NB; ++k) { b[k].key = 100u + (uint32_t)k; b[k].color.mode = MD_RENDER_COLOR_UNIFORM; }
    b[0].kind = MD_RENDER_SPHERES; b[0].count = 1; b[0].points.count = 1; b[0].points.position = f_buf(&fx, opq, 12);
    b[0].style.spheres.radius_add = 5.0f; b[0].color.uniform = 0xFF808080u;
    md_gpu_addr_t GA = f_buf(&fx, ga, sizeof(ga));
    b[1].kind = MD_RENDER_SPHERES; b[1].count = 2; b[1].points.count = 2; b[1].points.position = GA;
    b[1].style.spheres.radius_add = 0.9f; b[1].color.uniform = 0xFF2040E0u;
    b[2].kind = MD_RENDER_CAPSULES; b[2].count = 1; b[2].points.count = 2; b[2].points.position = GA;
    b[2].style.capsules.radius = 0.35f; b[2].color.uniform = 0xFF20E040u;
    for (int k = 1; k <= 2; ++k) { b[k].color.transparency = 0.4f; b[k].color.group = 77; }
    b[3].kind = MD_RENDER_TRIANGLES; b[3].count = s.nt; b[3].points.count = s.nv;
    b[3].points.position = f_buf(&fx, mpos, s.nv * 12); b[3].index = f_buf(&fx, s.tri, s.nt * 12);
    b[3].style.triangles.normal = f_buf(&fx, mnrm, s.nv * 12);
    b[3].color.uniform = 0xFF20A0E0u; b[3].color.transparency = 0.5f;
    static const uint32_t stc[6] = { 0xFF0000F0u, 0xFF00F000u, 0xFFF00000u, 0xFF00F0F0u, 0xFFF000F0u, 0xFFF0F000u };
    for (int j = 0; j < 6; ++j) {
        md_render_batch_t* x = &b[4 + j];
        x->kind = MD_RENDER_SPHERES; x->count = 1; x->points.count = 1; x->points.position = f_buf(&fx, st + 3 * j, 12);
        x->style.spheres.radius_add = 0.5f; x->color.uniform = stc[j]; x->color.transparency = 0.3f;
    }

    /* The reference */
    int nq = 0;
    tp_ref_t* q = calloc(16 + s.nt, sizeof(tp_ref_t));
    for (int j = 0; j < 2; ++j) { tp_ref_t* x = &q[nq++]; x->kind = 0; for (int i = 0; i < 3; ++i) x->a[i] = ga[3 * j + i]; x->r = 0.9; x->group = 77; x->color = b[1].color.uniform; x->opacity = 0.6; }
    { tp_ref_t* x = &q[nq++]; x->kind = 1; for (int i = 0; i < 3; ++i) { x->a[i] = ga[i]; x->b[i] = ga[3 + i]; } x->r = 0.35; x->group = 77; x->color = b[2].color.uniform; x->opacity = 0.6; }
    for (uint32_t e = 0; e < s.nt; ++e) {
        tp_ref_t* x = &q[nq++];
        x->kind = 2; x->group = 100000u + e; x->color = b[3].color.uniform; x->opacity = 0.5;     /* no group: every fragment */
        for (int j = 0; j < 3; ++j) for (int i = 0; i < 3; ++i) { x->tri.q[j][i] = mpos[3 * s.tri[3 * e + j] + i]; x->tri.n[j][i] = mnrm[3 * s.tri[3 * e + j] + i]; }
    }
    for (int j = 0; j < 6; ++j) { tp_ref_t* x = &q[nq++]; x->kind = 0; for (int i = 0; i < 3; ++i) x->a[i] = st[3 * j + i]; x->r = 0.5; x->group = 1004 + (uint32_t)j; x->color = stc[j]; x->opacity = 0.7; }
    tp_ref_t op = { .kind = 0, .a = { opq[0], opq[1], opq[2] }, .r = 5.0 };

    md_render_frame_t f = frame_of(W, H, cam, b, NB);
    md_render_status_t stt;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &stt));
    md_render_targets_t tg = md_render_view_targets(fx.view);
    ASSERT_TRUE(tg.layers != 0);
    ASSERT_EQ(K, tg.layer_count);
    ASSERT_EQ(W * H, tg.layer_stride);
    const size_t P = (size_t)W * H;
    uint32_t* L  = malloc(P * K * 12);
    uint32_t* id = malloc(P * 8);
    ASSERT_TRUE(f_read_buf(&fx, tg.layers, L, P * K * 12) && f_read(&fx, tg.id, id, P * 8));

    int bad = 0, sure = 0, deduped = 0, full = 0, ghost = 0, mesh_backs = 0, opaque_only = 0, nerr = 0;
    double worst_t = 0;
    double* errs = malloc(P * K * sizeof(double));
    for (uint32_t y = 0; y < H; ++y) {
        for (uint32_t x = 0; x < W; ++x) {
            const uint32_t p = y * W + x;
            /* the opaque ids: the background sphere's, or nothing; no transparent one */
            if (id[2 * p + 1] != 0 && (id[2 * p + 1] >> 20) != 1u) { bad++; if (bad < 8) printf("    pixel (%u, %u): a transparent surface in the ids\n", x, y); }
            int lay[8], nfr;
            double lt[8];
            ray_t r = pixel_ray(&c, W, H, x + 0.5, y + 0.5);
            double t_op, nop[3];
            if (!tp_hit(&op, r, &t_op, nop)) t_op = 1e30;
            int nl = tp_layers(q, nq, r, t_op, K, lay, lt, &nfr);
            if (nl < 0) continue;
            bool agree = true;
            static const double off[4][2] = { { 0.05, 0 }, { -0.05, 0 }, { 0, 0.05 }, { 0, -0.05 } };
            for (int o = 0; o < 4 && agree; ++o) {
                int l2[8], nf2;
                double t2[8], to2, n2[3];
                ray_t r2 = pixel_ray(&c, W, H, x + 0.5 + off[o][0], y + 0.5 + off[o][1]);
                if (!tp_hit(&op, r2, &to2, n2)) to2 = 1e30;
                int nl2 = tp_layers(q, nq, r2, to2, K, l2, t2, &nf2);
                agree = nl2 == nl && nf2 == nfr && (t_op < 1e29) == (to2 < 1e29);
                for (int j = 0; j < nl && agree; ++j) agree = q[l2[j]].group == q[lay[j]].group;
            }
            if (!agree) continue;
            sure++;
            deduped += nfr > nl && nl < (int)K;
            full    += nfr > (int)K;
            for (int j = 0; j < nl; ++j) ghost += q[lay[j]].group == 77;
            if (nl == 0 && t_op < 1e29) opaque_only++;
            for (int j = 0; j < (int)K; ++j) {
                const float  dl = *(const float*)&L[2 * K * P + j * P + p];
                if (j >= nl) {
                    if (j == nl && dl != 0.0f) { bad++; if (bad < 8) printf("    pixel (%u, %u): layer %d of %d has depth %g\n", x, y, j, nl, dl); }
                    continue;
                }
                const tp_ref_t* e = &q[lay[j]];
                /* depth: reverse Z, n / -z */
                const double tz = (c.n / dl) / -r.d[2];
                const double et = fabs(tz - lt[j]) / lt[j];
                worst_t = fmax(worst_t, et);
                const uint32_t col = L[j * P + p];
                const int ea = (int)lround(e->opacity * 255.0);
                int ce = 0;
                for (int ch = 0; ch < 3; ++ch) ce = MD_MAX(ce, abs((int)((col >> (8 * ch)) & 0xFF) - (int)((e->color >> (8 * ch)) & 0xFF)));
                if (et > 1e-4 || ce > 1 || abs((int)(col >> 24) - ea) > 1) {
                    bad++;
                    if (bad < 8) printf("    pixel (%u, %u) layer %d: depth %.5f against %.5f, colour %08x against %08x (alpha %d)\n", x, y, j, tz, lt[j], col, e->color, ea);
                    continue;
                }
                /* the normal, facing the eye */
                double t, n[3], ng[3];
                tp_hit(e, r, &t, n);
                if (dot3(n, r.d) > 0) { for (int i = 0; i < 3; ++i) n[i] = -n[i]; mesh_backs += e->kind == 2; }
                oct_decode(L[K * P + j * P + p], ng);
                double cs = dot3(n, ng);
                errs[nerr++] = acos(cs > 1 ? 1 : cs < -1 ? -1 : cs) * 180.0 / 3.14159265358979;
            }
        }
    }
    qsort(errs, (size_t)nerr, sizeof(double), cmp_double);
    const double p50 = nerr ? errs[nerr / 2] : 0, p99 = nerr ? errs[(nerr * 99) / 100] : 0;
    printf("    %d pixels sure: %d with fragments of a group dropped, %d with more than K, %d ghost layers, %d opaque alone; depth %.2e relative at most; normals over %d layers: median %.4f, 99th %.4f degrees\n",
           sure, deduped, full, ghost, opaque_only, worst_t, nerr, p50, p99);
    EXPECT_EQ(0, bad);
    EXPECT_GT(deduped, 50);
    EXPECT_GT(full, 3);
    EXPECT_GT(ghost, 200);
    EXPECT_GT(opaque_only, 200);
    EXPECT_LT(p50, 0.05);
    EXPECT_LT(p99, 1.0);
    (void)mesh_backs;

    /* The line of six: its pixel holds the four nearest, in order */
    {
        double pt[3] = { st[0], st[1], st[2] };
        const uint32_t px = (uint32_t)((pt[0] / -pt[2] / (c.tan_half * c.aspect) * 0.5 + 0.5) * W);
        const uint32_t py = (uint32_t)((0.5 - pt[1] / -pt[2] / c.tan_half * 0.5) * H);
        for (uint32_t j = 0; j < K; ++j) EXPECT_EQ(stc[j] & 0xFFFFFFu, L[j * P + py * W + px] & 0xFFFFFFu);
    }

    /* The two-pass fallback (32-bit atomics): the same layers, to the bit */
    {
        uint32_t* L2 = malloc(P * K * 12);
        md_render_debug_view(fx.view, MD_RENDER_DEBUG_KBUFFER_SPLIT, 0.0f);
        ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &stt));
        EXPECT_NE(0u, stt.work & MD_RENDER_WORK_RASTER);
        ASSERT_TRUE(f_read_buf(&fx, md_render_view_targets(fx.view).layers, L2, P * K * 12));
        size_t diff = 0;
        for (size_t i = 0; i < P * K * 3; ++i) {
            /* the colour and normal behind a layer of depth 0 are left as they were */
            const size_t j = i % (P * K), plane = i / (P * K);
            if (plane < 2 && *(const float*)&L[2 * P * K + j] == 0.0f) continue;
            diff += L[i] != L2[i];
        }
        printf("    two passes: %zu words of %zu differ\n", diff, P * K * 3);
        EXPECT_EQ((size_t)0, diff);
        md_render_debug_view(fx.view, 0, 0.0f);
        ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &stt));     /* the one pass again (the change rasterises) */
        free(L2);
    }

    /* Appearance alone: a new transparency is a resolve, and the layers'
       opacity follows */
    b[1].color.transparency = b[2].color.transparency = 0.75f;
    f = frame_of(W, H, cam, b, NB);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &stt));
    EXPECT_EQ(0u, stt.work & MD_RENDER_WORK_RASTER);
    EXPECT_NE(0u, stt.work & MD_RENDER_WORK_RESOLVE);
    ASSERT_TRUE(f_read_buf(&fx, tg.layers, L, P * K * 12));
    int quarter = 0;
    for (size_t p = 0; p < P; ++p) for (uint32_t j = 0; j < K; ++j) {
        uint32_t col = L[j * P + p];
        if (*(const float*)&L[2 * K * P + j * P + p] != 0.0f && (col & 0xFFFFFFu) == (b[1].color.uniform & 0xFFFFFFu)) quarter += (col >> 24) == 64u;
    }
    EXPECT_GT(quarter, 100);

    /* Opaque again: a raster, and no layers */
    for (int k = 1; k < NB; ++k) b[k].color.transparency = 0.0f;
    f = frame_of(W, H, cam, b, NB);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &stt));
    EXPECT_NE(0u, stt.work & MD_RENDER_WORK_RASTER);
    EXPECT_EQ(0u, (unsigned)md_render_view_targets(fx.view).layers);

    /* Refused: a transparent spline or polygon batch, more than 8 layers */
    {
        md_render_batch_t x = b[1];
        x.kind = MD_RENDER_POLYGONS; x.style.polygons.corners = 3; x.style.polygons.radius = 0.1f; x.count = 0;
        x.color.transparency = 0.5f;
        x.count = 1; x.points.count = 3;
        md_render_frame_t fr = frame_of(W, H, cam, &x, 1);
        EXPECT_FALSE(md_render_view_draw(fx.gfx, fx.view, &fr, NULL));
        fr = frame_of(W, H, cam, b, NB);
        fr.layers = 9;
        EXPECT_FALSE(md_render_view_draw(fx.gfx, fx.view, &fr, NULL));
    }
    free(errs); free(L); free(id); free(q); free(mpos); free(mnrm); free(s.xyz); free(s.tri);
    f_close(&fx);
}

/* ---- Compose ------------------------------------------------------------------------ */

/* md_compose.slang's shading, in C: the reference */
typedef struct { float x, y, z; } cf3;
static cf3 cf(float x, float y, float z) { cf3 r = { x, y, z }; return r; }
static cf3 cadd(cf3 a, cf3 b) { return cf(a.x + b.x, a.y + b.y, a.z + b.z); }
static cf3 cmul(cf3 a, float s) { return cf(a.x * s, a.y * s, a.z * s); }
static cf3 cmulv(cf3 a, cf3 b) { return cf(a.x * b.x, a.y * b.y, a.z * b.z); }
static float cdot(cf3 a, cf3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static cf3 cnorm(cf3 a) { return cmul(a, 1.0f / sqrtf(cdot(a, a))); }
static float csat(float x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
static float cmix(float a, float b, float t) { return a + (b - a) * t; }
static cf3 cmix3(cf3 a, cf3 b, float t) { return cf(cmix(a.x, b.x, t), cmix(a.y, b.y, t), cmix(a.z, b.z, t)); }
static float c_s2l(float c) { c = csat(c); return c < 0.04045f ? c / 12.92f : powf((c + 0.055f) / 1.055f, 2.4f); }
static cf3 c_srgb(cf3 c) { return cf(c_s2l(c.x), c_s2l(c.y), c_s2l(c.z)); }

/* What a surface reflects diffusely and specularly, and Fe (md_compose.slang's Lit) */
typedef struct { cf3 d, s; float Fe; } c_lit_t;

static c_lit_t c_shade(const md_compose_shading_t* sh, cf3 color, cf3 V, cf3 N, float ao) {
    c_lit_t r = { cf(0, 0, 0), cf(0, 0, 0), 0.0f };
    cf3 L = cnorm(cf(sh->light_dir[0], sh->light_dir[1], sh->light_dir[2]));
    cf3 H = cnorm(cadd(L, V));
    if (sh->mode == MD_COMPOSE_SHADING_BASIC) {
        float w = 0.5f + 0.5f * cdot(N, L), amb = sh->ambient, light = sh->light_strength;
        cf3 c = cmul(cf(csat(color.x), csat(color.y), csat(color.z)), amb + light * (1.0f - amb) * w * w);
        float m = fmaxf(1.0f, fmaxf(c.x, fmaxf(c.y, c.z)));
        c = cmul(c, 1.0f / m);
        float sp = light * sh->specular * powf(fmaxf(cdot(N, H), 0.0f), 48.0f);
        c = cmix3(c, cf(1, 1, 1), csat(sp));
        r.d = cmul(sh->linear_colors ? c : c_srgb(c), ao);
        return r;
    }
    const float PI = 3.1415926535f;
    float HdotV = csat(cdot(H, V)), NdotH = csat(cdot(N, H)), NdotL = csat(cdot(N, L)), NdotV = csat(cdot(N, V));
    float F0 = csat(sh->f0), rough = sh->roughness < 0.04f ? 0.04f : sh->roughness > 1 ? 1 : sh->roughness;
    cf3 albedo = sh->linear_colors ? cf(csat(color.x), csat(color.y), csat(color.z)) : c_srgb(color);
    cf3 sky = cf(sh->sky_radiance[0], sh->sky_radiance[1], sh->sky_radiance[2]);
    cf3 gnd = cf(sh->ground_radiance[0], sh->ground_radiance[1], sh->ground_radiance[2]);
    if (NdotL > 0) {
        float a = rough * rough, a2 = a * a, den = NdotH * NdotH * (a2 - 1) + 1, NDF = a2 / (PI * den * den);
        float rr = rough + 1, k = rr * rr / 8;
        float G = (NdotL / (NdotL * (1 - k) + k)) * (NdotV / (NdotV * (1 - k) + k));
        float F = F0 + (1 - F0) * powf(csat(1 - HdotV), 5.0f);
        float spec = NDF * G * F / (4 * NdotV * NdotL + 0.0001f);
        cf3 E = cmul(cf(sh->dir_radiance[0], sh->dir_radiance[1], sh->dir_radiance[2]), sh->light_strength * NdotL);
        r.s = cadd(r.s, cmul(E, spec));
        r.d = cadd(r.d, cmulv(cmul(albedo, (1 - F) / PI), E));
    }
    float rx = rough * -1.0f + 1.0f, ry = rough * -0.0275f + 0.0425f, rz = rough * -0.572f + 1.04f, rw = rough * 0.022f - 0.04f;
    float a004 = fminf(rx * rx, exp2f(-9.28f * NdotV)) * rx + ry;
    float abx = -1.04f * a004 + rz, aby = 1.04f * a004 + rw;
    r.Fe = F0 * abx + aby;
    cf3 R = cadd(cmul(V, -1), cmul(N, 2 * cdot(V, N)));         /* reflect(-V, N) */
    float wv = 0.1f + rough, tt = csat((R.y + wv) / (2 * wv)), smooth = tt * tt * (3 - 2 * tt);
    const float so = csat(powf(NdotV + ao, exp2f(-16.0f * rough * rough - 1.0f)) - 1.0f + ao);    /* specular occlusion */
    r.d = cadd(r.d, cmul(cmulv(albedo, cmix3(gnd, sky, 0.5f + 0.5f * N.y)), (1 - r.Fe) * ao));
    r.s = cadd(r.s, cmul(cmix3(gnd, sky, smooth), r.Fe * so));
    return r;
}

/* The scene of a view, both models, without occlusion: the opaque surface
   and the layers over it against the reference over the same targets and
   layers (read back) */
UTEST(render, compose_shades_the_scene) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    md_gpu_device_info_t info;
    md_gpu_device_info(fx.dev, &info);
    if (!info.fragment_stores) UTEST_SKIP("fragment stages cannot write memory (the K-buffer's)");
    const uint32_t W = 96, H = 72, K = 4;
    cam_t c = { .tan_half = 0.4f, .aspect = (float)W / H, .n = 0.1f };
    md_render_camera_t cam = make_camera(&c);
    ico_t s = ico_make(2);
    float* mpos = malloc(s.nv * 12);
    for (uint32_t i = 0; i < s.nv; ++i) for (int k = 0; k < 3; ++k) mpos[3 * i + k] = (k == 2 ? -13.0f : k == 0 ? 0.6f : -0.3f) + 1.7f * s.xyz[3 * i + k];
    float opq[3] = { 0.0f, 0.0f, -22.0f };
    float ga[6] = { -2.0f, 0.5f, -12.0f,   0.0f, 1.5f, -12.5f };
    md_render_batch_t b[4];
    memset(b, 0, sizeof(b));
    for (int k = 0; k < 4; ++k) { b[k].key = 200u + (uint32_t)k; b[k].color.mode = MD_RENDER_COLOR_UNIFORM; }
    b[0].kind = MD_RENDER_SPHERES; b[0].count = 1; b[0].points.count = 1; b[0].points.position = f_buf(&fx, opq, 12);
    b[0].style.spheres.radius_add = 5.0f; b[0].color.uniform = 0xFF8060A0u;
    md_gpu_addr_t GA = f_buf(&fx, ga, sizeof(ga));
    b[1].kind = MD_RENDER_SPHERES; b[1].count = 2; b[1].points.count = 2; b[1].points.position = GA;
    b[1].style.spheres.radius_add = 0.9f; b[1].color.uniform = 0xFF2040E0u;
    b[2].kind = MD_RENDER_CAPSULES; b[2].count = 1; b[2].points.count = 2; b[2].points.position = GA;
    b[2].style.capsules.radius = 0.35f; b[2].color.uniform = 0xFF20E040u;
    for (int k = 1; k <= 2; ++k) { b[k].color.transparency = 0.4f; b[k].color.group = 77; }
    b[3].kind = MD_RENDER_TRIANGLES; b[3].count = s.nt; b[3].points.count = s.nv;
    b[3].points.position = f_buf(&fx, mpos, s.nv * 12); b[3].index = f_buf(&fx, s.tri, s.nt * 12);
    b[3].color.uniform = 0xFF20A0E0u; b[3].color.transparency = 0.5f;
    md_render_frame_t f = frame_of(W, H, cam, b, 4);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    md_render_targets_t tg = md_render_view_targets(fx.view);
    const size_t P = (size_t)W * H;
    uint32_t* L  = malloc(P * K * 12);
    uint32_t* oc = malloc(P * 4);
    uint32_t* on = malloc(P * 4);
    ASSERT_TRUE(f_read_buf(&fx, tg.layers, L, P * K * 12) && f_read(&fx, tg.color, oc, P * 4) && f_read(&fx, tg.normal, on, P * 4));

    md_compose_ctx_t cx = md_compose_create(fx.dev);
    ASSERT_TRUE(cx != NULL);
    md_gpu_texture_desc_t td = {0};
    td.type = MD_GPU_TEX_2D; td.usage = MD_GPU_TEX_STORAGE; td.width = W; td.height = H;
    td.format = MD_GPU_FORMAT_RGBA16_FLOAT; td.label = "compose";
    md_gpu_texture_t out = md_gpu_texture_create(fx.gfx, &td);
    uint16_t* px = malloc(P * 8);
    for (int mode = 0; mode < 2; ++mode) {
        md_compose_desc_t d = {0};
        d.camera = cam;
        d.target = out;
        md_compose_shading_t* sh = &d.shading;
        sh->mode = mode ? MD_COMPOSE_SHADING_BASIC : MD_COMPOSE_SHADING_PBR;
        sh->linear_colors = mode == 1;
        sh->light_dir[0] = 0.3f; sh->light_dir[1] = 0.8f; sh->light_dir[2] = 0.5f;
        sh->light_strength = 1.2f;
        for (int i = 0; i < 3; ++i) { sh->dir_radiance[i] = 2.5f; sh->sky_radiance[i] = 0.4f + 0.1f * i; sh->ground_radiance[i] = 0.15f; }
        sh->roughness = 0.4f; sh->f0 = 0.04f; sh->ambient = 0.3f; sh->specular = 0.25f;
        ASSERT_TRUE(md_compose(fx.gfx, cx, fx.view, &d));
        ASSERT_TRUE(f_read(&fx, out, px, P * 8));
        int bad = 0, covered = 0, multi = 0, opaque = 0, both = 0;
        double worst = 0;
        for (uint32_t y = 0; y < H; ++y) {
            for (uint32_t x = 0; x < W; ++x) {
                const size_t p = (size_t)y * W + x;
                ray_t r = pixel_ray(&c, W, H, x + 0.5, y + 0.5);
                cf3 V = cf((float)-r.d[0], (float)-r.d[1], (float)-r.d[2]);
                cf3 C = cf(0, 0, 0);
                float T = 1.0f;
                int n = 0;
                for (uint32_t k = 0; k < K; ++k) {
                    if (*(const float*)&L[2 * K * P + k * P + p] == 0.0f) break;
                    const uint32_t c8 = L[k * P + p];
                    cf3 col = cf((c8 & 0xFF) / 255.0f, ((c8 >> 8) & 0xFF) / 255.0f, ((c8 >> 16) & 0xFF) / 255.0f);
                    const float a = (c8 >> 24) / 255.0f;
                    double nd[3];
                    oct_decode(L[K * P + k * P + p], nd);
                    cf3 N = cnorm(cf((float)nd[0], (float)nd[1], (float)nd[2]));
                    c_lit_t l = c_shade(sh, col, V, N, 1.0f);
                    C = cadd(C, cmul(cadd(cmul(l.d, a), l.s), T));        /* a film: a x diffuse + specular */
                    T *= (1.0f - a) * (1.0f - l.Fe);
                    ++n;
                }
                /* the opaque surface behind them */
                const float cov = (oc[p] >> 24) ? 1.0f : 0.0f;
                if (cov > 0) {
                    cf3 col = cf((oc[p] & 0xFF) / 255.0f, ((oc[p] >> 8) & 0xFF) / 255.0f, ((oc[p] >> 16) & 0xFF) / 255.0f);
                    double nd[3];
                    oct_decode(on[p], nd);
                    cf3 N = cnorm(cf((float)nd[0], (float)nd[1], (float)nd[2]));
                    c_lit_t l = c_shade(sh, col, V, N, 1.0f);
                    C = cadd(C, cmul(cadd(l.d, l.s), T));
                }
                covered += n > 0;
                multi   += n > 1;
                opaque  += cov > 0;
                both    += n > 0 && cov > 0;
                const float got[4] = { half_bits_to_float(px[4 * p]), half_bits_to_float(px[4 * p + 1]), half_bits_to_float(px[4 * p + 2]), half_bits_to_float(px[4 * p + 3]) };
                const float exp[4] = { C.x, C.y, C.z, 1.0f - T * (1.0f - cov) };
                for (int ch = 0; ch < 4; ++ch) {
                    const double e = fabs(got[ch] - exp[ch]) / (2e-3 * fmax(1.0, fabs(exp[ch])) + 2e-3);
                    worst = fmax(worst, e);
                    if (e > 1.0) { bad++; if (bad < 6) printf("    %s pixel (%u, %u) channel %d: %g against %g\n", mode ? "basic" : "pbr", x, y, ch, got[ch], exp[ch]); }
                }
            }
        }
        printf("    %s: %d pixels opaque, %d with layers (%d over the opaque surface), %d with more than one; at most %.2f of the tolerance\n",
               mode ? "basic" : "pbr", opaque, covered, both, multi, worst);
        EXPECT_EQ(0, bad);
        EXPECT_GT(covered, 300);
        EXPECT_GT(multi, 100);
        EXPECT_GT(opaque, 500);
        EXPECT_GT(both, 100);
    }
    md_gpu_texture_destroy(out);
    md_compose_destroy(cx);
    free(px); free(L); free(oc); free(on); free(mpos); free(s.xyz); free(s.tri);
    f_close(&fx);
}

/* The occlusion, by what it must do rather than how (the estimator's own
   numbers are viamd's): white surfaces under BASIC shading with no light and
   an ambient of 1 show the visibility itself. A sphere resting on a plane,
   smaller on the screen than the occlusion's reach (a tenth of the height):
   the plane is open far from it and occluded around the contact; the
   sphere's face toward the eye is open. Scale free: the scene and the camera
   scaled by 40 give the same visibility. Perspective and orthographic. An
   intensity of 0 estimates nothing (1 everywhere). */
static bool occ_draw(fix_t* fx, md_render_camera_t cam, uint32_t W, uint32_t H, float k, md_compose_ctx_t cx, md_gpu_texture_t out, float ao, md_compose_ao_mode_t mode, uint16_t* px) {
    const float z0 = -20.0f;
    float quad[12] = { -30, -30, z0,  30, -30, z0,  30, 30, z0,  -30, 30, z0 };
    const uint32_t tri[6] = { 0, 1, 2, 0, 2, 3 };
    float ball[3] = { 0.0f, 0.0f, z0 + 0.5f };
    for (int i = 0; i < 12; ++i) quad[i] *= k;
    for (int i = 0; i < 3; ++i) ball[i] *= k;
    md_render_batch_t b[2];
    memset(b, 0, sizeof(b));
    b[0].key = 1; b[0].kind = MD_RENDER_TRIANGLES; b[0].count = 2; b[0].points.count = 4;
    b[0].points.position = f_buf(fx, quad, sizeof(quad)); b[0].index = f_buf(fx, tri, sizeof(tri));
    b[1].key = 2; b[1].kind = MD_RENDER_SPHERES; b[1].count = 1; b[1].points.count = 1;
    b[1].points.position = f_buf(fx, ball, sizeof(ball)); b[1].style.spheres.radius_add = 0.5f * k;
    for (int i = 0; i < 2; ++i) { b[i].color.mode = MD_RENDER_COLOR_UNIFORM; b[i].color.uniform = 0xFFFFFFFFu; }
    md_render_frame_t f = frame_of(W, H, cam, b, 2);
    f.jitter = MD_RENDER_JITTER_NONE;
    if (!md_render_view_draw(fx->gfx, fx->view, &f, NULL)) return false;
    md_compose_desc_t d = {0};
    d.camera = cam;
    d.target = out;
    d.ao = ao;
    d.ao_mode = mode;
    d.shading.mode = MD_COMPOSE_SHADING_BASIC;
    d.shading.linear_colors = true;
    d.shading.light_dir[2] = 1.0f;
    d.shading.ambient = 1.0f;
    return md_compose(fx->gfx, cx, fx->view, &d) && f_read(fx, out, px, (size_t)W * H * 8);
}

UTEST(render, compose_occlusion) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 256, H = 192;
    const size_t P = (size_t)W * H;
    md_compose_ctx_t cx = md_compose_create(fx.dev);
    ASSERT_TRUE(cx != NULL);
    md_gpu_texture_desc_t td = {0};
    td.type = MD_GPU_TEX_2D; td.usage = MD_GPU_TEX_STORAGE; td.width = W; td.height = H;
    td.format = MD_GPU_FORMAT_RGBA16_FLOAT; td.label = "compose";
    md_gpu_texture_t out = md_gpu_texture_create(fx.gfx, &td);
    uint16_t* v1 = malloc(P * 8);
    uint16_t* v2 = malloc(P * 8);
    for (int run = 0; run < 4; ++run) {
        const int ortho = run & 1;
        const md_compose_ao_mode_t mode = run < 2 ? MD_COMPOSE_AO_QUALITY : MD_COMPOSE_AO_PERFORMANCE;
        const char* mname = mode == MD_COMPOSE_AO_QUALITY ? "quality" : "performance";
        const float k = 40.0f;
        cam_t c1 = ortho ? (cam_t){ .ortho = true, .half_h = 6.0f, .aspect = (float)W / H, .n = 0.0f, .far_ = 100.0f }
                         : (cam_t){ .tan_half = 0.3f, .aspect = (float)W / H, .n = 0.1f };
        cam_t c2 = c1;
        if (ortho) { c2.half_h *= k; c2.far_ *= k; } else { c2.n *= k; }
        ASSERT_TRUE(occ_draw(&fx, make_camera(&c1), W, H, 1.0f, cx, out, 5.0f, mode, v1));
        ASSERT_TRUE(occ_draw(&fx, make_camera(&c2), W, H, k, cx, out, 5.0f, mode, v2));
        /* the sphere's silhouette: its radius on the screen, about its centre */
        const double rs = ortho ? 0.5 / c1.half_h * 0.5 * H : 0.5 / (19.5 * c1.tan_half) * 0.5 * H;
        const double cx0 = 0.5 * W, cy0 = 0.5 * H;
        double far_min = 1, ring_min = 1, ring_sum = 0, top_min = 1, worst = 0;
        int ring_n = 0, far_n = 0, bad = 0;
        for (uint32_t y = 0; y < H; ++y) {
            for (uint32_t x = 0; x < W; ++x) {
                const size_t p = (size_t)y * W + x;
                const double v = half_bits_to_float(v1[4 * p]), w = half_bits_to_float(v2[4 * p]);
                if (!(v >= 0.0 && v <= 1.0)) { bad++; continue; }
                worst = fmax(worst, fabs(v - w));
                const double r = sqrt((x + 0.5 - cx0) * (x + 0.5 - cx0) + (y + 0.5 - cy0) * (y + 0.5 - cy0));
                if (r > rs + 0.1 * H + 16.0) { far_min = fmin(far_min, v); far_n++; }   /* beyond its reach and the blur's */
                else if (r > rs + 1.5 && r < rs + 4.0) { ring_min = fmin(ring_min, v); ring_sum += v; ring_n++; }
                else if (r < 0.3 * rs) top_min = fmin(top_min, v);
            }
        }
        printf("    %s, %s: far from the sphere %.3f at least (%d px), around its contact %.3f mean, %.3f least (%d px); its face %.3f at least; scaled by %g: %.4f off at most\n",
               mname, ortho ? "orthographic" : "perspective", far_min, far_n, ring_sum / fmax(ring_n, 1), ring_min, ring_n, top_min, k, worst);
        EXPECT_EQ(0, bad);
        EXPECT_GT(far_min, 0.97);
        EXPECT_GT(top_min, 0.97);
        EXPECT_LT(ring_sum / fmax(ring_n, 1), 0.92);
        EXPECT_LT(worst, 0.02);
        /* the intensity: twice of it, twice the obscurance (performance, exactly: the blur is linear in
           it) or the visibility squared (quality, before the blur: after it, within its spread) */
        ASSERT_TRUE(occ_draw(&fx, make_camera(&c1), W, H, 1.0f, cx, out, 10.0f, mode, v2));
        double lin = 0;
        for (size_t p = 0; p < P; ++p) {
            const double v = half_bits_to_float(v1[4 * p]), w = half_bits_to_float(v2[4 * p]);
            lin = fmax(lin, fabs(w - (mode == MD_COMPOSE_AO_PERFORMANCE ? 1.0 - 2.0 * (1.0 - v) : v * v)));
        }
        printf("    intensity doubled: %.4f off at most\n", lin);
        EXPECT_LT(lin, mode == MD_COMPOSE_AO_PERFORMANCE ? 4e-3 : 0.03);
        EXPECT_GT(far_n, 1000);
        EXPECT_GT(ring_n, 50);
    }
    /* no occlusion asked: none */
    {
        cam_t c1 = { .tan_half = 0.3f, .aspect = (float)W / H, .n = 0.1f };
        ASSERT_TRUE(occ_draw(&fx, make_camera(&c1), W, H, 1.0f, cx, out, 0.0f, MD_COMPOSE_AO_QUALITY, v1));
        int open = 0, covered = 0;
        for (size_t p = 0; p < P; ++p) {
            if (half_bits_to_float(v1[4 * p + 3]) == 0.0f) continue;
            covered++;
            open += half_bits_to_float(v1[4 * p]) == 1.0f;
        }
        EXPECT_EQ(covered, open);
        EXPECT_GT(covered, (int)(P / 2));
    }
    md_gpu_texture_destroy(out);
    md_compose_destroy(cx);
    free(v1); free(v2);
    f_close(&fx);
}

static uint16_t float_to_half(float f);

/* ---- Volumes in the compose ---------------------------------------------------------- */

static md_gpu_texture_t f_tex(fix_t* fx, md_gpu_format_t format, md_gpu_tex_usage_t usage, uint32_t w, uint32_t h, const char* label);
static bool f_upload_tex(fix_t* fx, md_gpu_texture_t t, const void* src, size_t size);

/* A box [-2, 2]^2 x [-14, -10] of world: texture space [0, 1]^3 scaled by 4 */
static void vol_model(float m[16]) {
    memset(m, 0, 16 * sizeof(float));
    m[0] = m[5] = m[10] = 4.0f;
    m[12] = -2.0f; m[13] = -2.0f; m[14] = -14.0f; m[15] = 1.0f;
}

static md_gpu_texture_t f_tex3d(fix_t* fx, uint32_t n, const float* v) {
    md_gpu_texture_desc_t td = {0};
    td.type = MD_GPU_TEX_3D; td.usage = MD_GPU_TEX_SAMPLED; td.width = td.height = td.depth_or_layers = n;
    td.format = MD_GPU_FORMAT_R32_FLOAT; td.label = "volume";
    md_gpu_texture_t t = md_gpu_texture_create(fx->gfx, &td);
    md_gpu_tex_region_t r = {0};
    if (t && !md_gpu_upload_texture(fx->gfx, t, &r, v, (size_t)n * n * n * 4)) return 0;
    return t;
}

/* The interval of the ray inside the sphere (c, r); false if it misses */
static bool ray_sphere(ray_t r, const double c[3], double rad, double* t0, double* t1) {
    double oc[3] = { r.o[0] - c[0], r.o[1] - c[1], r.o[2] - c[2] };
    double b = dot3(oc, r.d), cc = dot3(oc, oc) - rad * rad, disc = b * b - cc;
    if (disc < 0) return false;
    *t0 = -b - sqrt(disc); *t1 = -b + sqrt(disc);
    return true;
}

/* The interval of the ray inside the box [lo, hi] */
static bool ray_box(ray_t r, const double lo[3], const double hi[3], double* t0, double* t1) {
    double a = -1e30, b = 1e30;
    for (int k = 0; k < 3; ++k) {
        if (fabs(r.d[k]) < 1e-15) { if (r.o[k] < lo[k] || r.o[k] > hi[k]) return false; continue; }
        double u = (lo[k] - r.o[k]) / r.d[k], v = (hi[k] - r.o[k]) / r.d[k];
        a = fmax(a, fmin(u, v)); b = fmin(b, fmax(u, v));
    }
    *t0 = fmax(a, 0.0); *t1 = b;
    return *t0 < *t1;
}

static void unlit_shading(md_compose_shading_t* sh) {   /* BASIC, no light, ambient 1, linear: a colour shows as it is */
    memset(sh, 0, sizeof(*sh));
    sh->mode = MD_COMPOSE_SHADING_BASIC;
    sh->linear_colors = true;
    sh->light_dir[2] = 1.0f;
    sh->ambient = 1.0f;
}

/* Volumes composited with what md_render drew, against analytic references:
   - emission of a constant field through a constant transfer function: the
     radiance c (1 - T) and the transmittance (1 - alpha)^(150 x the path in
     texture space), cut by the box, a plane and an opaque sphere inside it
     (white, seen through what is left), perspective and orthographic
   - an isosurface of a distance field: a sphere, its shading against the
     analytic sphere's at the analytic hit
   - optical density: an invisible isosurface (coverage 0) absorbing outside a
     ball, around an opaque sphere inside the ball, then around a transparent
     one (its front a layer met inside the volume) */
UTEST(render, compose_volumes) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    md_gpu_device_info_t info;
    md_gpu_device_info(fx.dev, &info);
    const bool layers_ok = info.fragment_stores;
    const uint32_t W = 96, H = 72, P = W * H, N = 48;
    const double box_lo[3] = { -2, -2, -14 }, box_hi[3] = { 2, 2, -10 }, ctr[3] = { 0, 0, -12 };
    md_compose_ctx_t cx = md_compose_create(fx.dev);
    ASSERT_TRUE(cx != NULL);
    md_gpu_texture_t out = f_tex(&fx, MD_GPU_FORMAT_RGBA16_FLOAT, MD_GPU_TEX_STORAGE, W, H, "compose");
    uint16_t* px = malloc(P * 8);

    /* the fields: constant, and the distance from the box's centre in world units */
    float* fc = malloc(N * N * N * 4);
    float* fd = malloc(N * N * N * 4);
    for (uint32_t z = 0; z < N; ++z) for (uint32_t y = 0; y < N; ++y) for (uint32_t x = 0; x < N; ++x) {
        const size_t i = ((size_t)z * N + y) * N + x;
        const double u[3] = { (x + 0.5) / N - 0.5, (y + 0.5) / N - 0.5, (z + 0.5) / N - 0.5 };
        fc[i] = 0.5f;
        fd[i] = (float)(4.0 * sqrt(dot3(u, u)));
    }
    md_gpu_texture_t vc = f_tex3d(&fx, N, fc), vdist = f_tex3d(&fx, N, fd);
    ASSERT_TRUE(vc && vdist);
    const uint8_t tfc[3] = { 200, 120, 40 }, tfa = 3;
    uint32_t tf[16];
    for (int i = 0; i < 16; ++i) tf[i] = tfc[0] | (tfc[1] << 8) | (tfc[2] << 16) | ((uint32_t)tfa << 24);
    md_gpu_texture_t tft = f_tex(&fx, MD_GPU_FORMAT_RGBA8_UNORM, MD_GPU_TEX_SAMPLED, 16, 1, "transfer");
    ASSERT_TRUE(f_upload_tex(&fx, tft, tf, sizeof(tf)));
    const cf3 tfl = c_srgb(cf(tfc[0] / 255.0f, tfc[1] / 255.0f, tfc[2] / 255.0f));
    const double alpha = tfa / 255.0;

    /* an opaque white sphere inside the box, or nothing */
    const float opq[3] = { 0.3f, -0.2f, -11.8f };
    const double opq_c[3] = { opq[0], opq[1], opq[2] }, opq_r = 0.7;
    md_render_batch_t ob = {0};
    ob.key = 1; ob.kind = MD_RENDER_SPHERES; ob.count = 1; ob.points.count = 1; ob.points.position = f_buf(&fx, opq, 12);
    ob.style.spheres.radius_add = (float)opq_r; ob.color.mode = MD_RENDER_COLOR_UNIFORM; ob.color.uniform = 0xFFFFFFFFu;

    for (int ortho = 0; ortho < 2; ++ortho) {
        cam_t c = ortho ? (cam_t){ .ortho = true, .half_h = 3.0f, .aspect = (float)W / H, .n = 0.0f, .far_ = 100.0f }
                        : (cam_t){ .tan_half = 0.25f, .aspect = (float)W / H, .n = 0.1f };
        md_render_camera_t cam = make_camera(&c);
        for (int cut = 0; cut < 2; ++cut) {
            md_render_frame_t f = frame_of(W, H, cam, &ob, 1);
            ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
            md_compose_volume_t v = {0};
            v.field = vc;
            vol_model(v.model);
            v.clip_max[0] = v.clip_max[1] = v.clip_max[2] = 1.0f;
            v.transfer = tft; v.transfer_range[0] = 0.0f; v.transfer_range[1] = 1.0f;
            if (cut) { v.plane_count = 1; v.planes[0][0] = 1.0f; v.planes[0][3] = -0.4f; }   /* x >= 0.4 kept */
            md_compose_desc_t d = {0};
            d.camera = cam; d.target = out; d.volumes = &v; d.volume_count = 1; d.frame = 3;
            unlit_shading(&d.shading);
            ASSERT_TRUE(md_compose(fx.gfx, cx, fx.view, &d));
            ASSERT_TRUE(f_read(&fx, out, px, P * 8));
            int bad = 0, through = 0, onto = 0;
            double worst = 0;
            for (uint32_t y = 0; y < H; ++y) for (uint32_t x = 0; x < W; ++x) {
                const size_t p = (size_t)y * W + x;
                ray_t r = pixel_ray(&c, W, H, x + 0.5, y + 0.5);
                double t0, t1, s0, s1;
                bool in = ray_box(r, box_lo, box_hi, &t0, &t1);
                const bool hit = ray_sphere(r, opq_c, opq_r, &s0, &s1) && s0 > 0;
                if (cut && in) {        /* x >= -2 + 0.4 x 4 = -0.4 */
                    if (fabs(r.d[0]) < 1e-12) { if (r.o[0] < -0.4) in = false; }
                    else { const double tx = (-0.4 - r.o[0]) / r.d[0]; if (r.d[0] > 0) t0 = fmax(t0, tx); else t1 = fmin(t1, tx); in = t0 < t1; }
                }
                if (in && hit) t1 = fmin(t1, s0);
                const double len = in ? fmax(t1 - t0, 0.0) / 4.0 : 0.0;       /* texture space */
                const double T = pow(1.0 - alpha, 150.0 * len);
                const double e[4] = { tfl.x * (1 - T) + (hit ? T : 0), tfl.y * (1 - T) + (hit ? T : 0), tfl.z * (1 - T) + (hit ? T : 0), hit ? 1.0 : 1 - T };
                /* not judged: a pixel the box's or the plane's or the sphere's edge crosses */
                bool edge = false;
                for (int k = 0; k < 4 && !edge; ++k) {
                    ray_t q = pixel_ray(&c, W, H, x + 0.5 + (k & 1 ? 0.6 : -0.6), y + 0.5 + (k & 2 ? 0.6 : -0.6));
                    double a0, a1, b0, b1;
                    bool qin = ray_box(q, box_lo, box_hi, &a0, &a1), qhit = ray_sphere(q, opq_c, opq_r, &b0, &b1);
                    if (cut && qin) qin = qin && (q.o[0] + q.d[0] * a1 >= -0.4 || q.o[0] + q.d[0] * a0 >= -0.4);
                    const bool qcut = cut && ((q.o[0] + q.d[0] * a0 < -0.4) != (r.o[0] + r.d[0] * t0 < -0.4));
                    edge = qin != in || qhit != hit || qcut;
                }
                if (edge) continue;
                through += in && !hit;
                onto    += in && hit;
                for (int k = 0; k < 4; ++k) {
                    const double err = fabs(half_bits_to_float(px[4 * p + k]) - e[k]);
                    worst = fmax(worst, err);
                    if (err > 4e-3) { bad++; if (bad < 6) printf("    emission%s%s pixel (%u, %u) channel %d: %g against %g\n", ortho ? ", orthographic" : "", cut ? ", cut" : "", x, y, k, half_bits_to_float(px[4 * p + k]), e[k]); }
                }
            }
            printf("    emission%s%s: %d pixels through the volume, %d onto the sphere in it; %.4f off at most\n",
                   ortho ? ", orthographic" : "", cut ? ", cut by a plane" : "", through, onto, worst);
            EXPECT_EQ(0, bad);
            EXPECT_GT(through, 300);
            EXPECT_GT(onto, 100);
        }

        /* an isosurface of the distance: a sphere of radius 1.3 about the box's centre, opaque, lit */
        {
            md_render_frame_t f = frame_of(W, H, cam, NULL, 0);
            ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
            md_compose_volume_t v = {0};
            v.field = vdist;
            vol_model(v.model);
            v.clip_max[0] = v.clip_max[1] = v.clip_max[2] = 1.0f;
            v.iso_count = 1; v.iso_values[0] = 1.3f;
            v.iso_colors[0][0] = 0.9f; v.iso_colors[0][1] = 0.5f; v.iso_colors[0][2] = 0.2f; v.iso_colors[0][3] = 1.0f;
            md_compose_desc_t d = {0};
            d.camera = cam; d.target = out; d.volumes = &v; d.volume_count = 1;
            md_compose_shading_t* sh = &d.shading;
            sh->mode = MD_COMPOSE_SHADING_PBR; sh->linear_colors = true;
            sh->light_dir[0] = 0.4f; sh->light_dir[1] = 0.6f; sh->light_dir[2] = 0.7f; sh->light_strength = 1.0f;
            for (int i = 0; i < 3; ++i) { sh->dir_radiance[i] = 2.7f; sh->sky_radiance[i] = 0.59f; sh->ground_radiance[i] = 0.49f; }
            sh->roughness = 0.4f; sh->f0 = 0.04f;
            ASSERT_TRUE(md_compose(fx.gfx, cx, fx.view, &d));
            ASSERT_TRUE(f_read(&fx, out, px, P * 8));
            int bad = 0, covered = 0;
            double worst = 0;
            for (uint32_t y = 0; y < H; ++y) for (uint32_t x = 0; x < W; ++x) {
                const size_t p = (size_t)y * W + x;
                ray_t r = pixel_ray(&c, W, H, x + 0.5, y + 0.5);
                double t0, t1, u0, u1;
                const bool hit = ray_sphere(r, ctr, 1.3, &t0, &t1);
                /* a voxel (1/12) from the silhouette either way */
                const bool inner = ray_sphere(r, ctr, 1.3 - 0.12, &u0, &u1), outer = ray_sphere(r, ctr, 1.3 + 0.12, &u0, &u1);
                if (inner != outer) continue;
                const double a = half_bits_to_float(px[4 * p + 3]);
                if (!hit) { if (a != 0.0) bad++; continue; }
                covered++;
                double hp[3], n[3];
                for (int k = 0; k < 3; ++k) hp[k] = r.o[k] + r.d[k] * t0;
                for (int k = 0; k < 3; ++k) n[k] = (hp[k] - ctr[k]) / 1.3;
                c_lit_t l = c_shade(sh, cf(0.9f, 0.5f, 0.2f), cf((float)-r.d[0], (float)-r.d[1], (float)-r.d[2]), cf((float)n[0], (float)n[1], (float)n[2]), 1.0f);
                const double e[4] = { l.d.x + l.s.x, l.d.y + l.s.y, l.d.z + l.s.z, 1.0 };
                for (int k = 0; k < 4; ++k) {
                    const double err = fabs(half_bits_to_float(px[4 * p + k]) - e[k]);
                    worst = fmax(worst, err);
                    if (err > 0.02) { bad++; if (bad < 6) printf("    isosurface%s pixel (%u, %u) channel %d: %g against %g\n", ortho ? ", orthographic" : "", x, y, k, half_bits_to_float(px[4 * p + k]), e[k]); }
                }
            }
            printf("    isosurface%s: %d pixels on it, against the analytic sphere %.4f off at most\n", ortho ? ", orthographic" : "", covered, worst);
            EXPECT_EQ(0, bad);
            EXPECT_GT(covered, 300);
        }

        /* optical density: an isosurface of coverage 0 (radius 1.5) absorbing,
           around an opaque white sphere, then around a transparent one */
        for (int inner_layers = 0; inner_layers < 2; ++inner_layers) {
            if (inner_layers && !layers_ok) break;
            const double sig = 0.6, R = 1.5, rs = 0.6, ta = 0.5;
            const float sp[3] = { (float)ctr[0], (float)ctr[1], (float)ctr[2] };
            md_render_batch_t b = {0};
            b.key = 7; b.kind = MD_RENDER_SPHERES; b.count = 1; b.points.count = 1; b.points.position = f_buf(&fx, sp, 12);
            b.style.spheres.radius_add = (float)rs; b.color.mode = MD_RENDER_COLOR_UNIFORM; b.color.uniform = 0xFFFFFFFFu;
            if (inner_layers) b.color.transparency = (float)(1.0 - ta);
            md_render_frame_t f = frame_of(W, H, cam, &b, 1);
            ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
            md_compose_volume_t v = {0};
            v.field = vdist;
            vol_model(v.model);
            v.clip_max[0] = v.clip_max[1] = v.clip_max[2] = 1.0f;
            v.iso_count = 1; v.iso_values[0] = (float)R; v.iso_densities[0] = (float)sig;
            v.iso_colors[0][0] = v.iso_colors[0][1] = v.iso_colors[0][2] = 1.0f;     /* coverage 0 */
            md_compose_desc_t d = {0};
            d.camera = cam; d.target = out; d.volumes = &v; d.volume_count = 1;
            unlit_shading(&d.shading);
            ASSERT_TRUE(md_compose(fx.gfx, cx, fx.view, &d));
            ASSERT_TRUE(f_read(&fx, out, px, P * 8));
            int bad = 0, judged = 0;
            double worst = 0;
            for (uint32_t y = 0; y < H; ++y) for (uint32_t x = 0; x < W; ++x) {
                const size_t p = (size_t)y * W + x;
                ray_t r = pixel_ray(&c, W, H, x + 0.5, y + 0.5);
                double a0, a1, b0, b1, u0, u1;
                /* "inside": the distance field's value at or above the isovalue is inside (v >= 0), so the
                   absorbing part is OUTSIDE radius R within the box; the ray's path in the box minus the ball */
                double t0, t1;
                if (!ray_box(r, box_lo, box_hi, &t0, &t1)) continue;
                bool box_edge = false;
                for (int k = 0; k < 4; ++k) {
                    ray_t q = pixel_ray(&c, W, H, x + 0.5 + (k & 1 ? 0.6 : -0.6), y + 0.5 + (k & 2 ? 0.6 : -0.6));
                    double q0, q1;
                    box_edge |= !ray_box(q, box_lo, box_hi, &q0, &q1);
                }
                if (box_edge) continue;
                const bool ball = ray_sphere(r, ctr, R, &a0, &a1);
                const bool hit  = ray_sphere(r, ctr, rs, &b0, &b1);
                /* not judged near the silhouettes */
                if (ray_sphere(r, ctr, R - 0.1, &u0, &u1) != ray_sphere(r, ctr, R + 0.1, &u0, &u1)) continue;
                if (ray_sphere(r, ctr, rs - 0.05, &u0, &u1) != ray_sphere(r, ctr, rs + 0.05, &u0, &u1)) continue;
                double C = 0, T = 1, t = t0;
                /* the absorbing stretches: [t0, a0] and [a1, t1] when the ray meets the ball, else [t0, t1] */
                #define ABSORB_TO(tt) do { double e0 = t, e1 = (tt); if (ball) { double s = fmax(0, fmin(e1, a0) - e0) + fmax(0, e1 - fmax(e0, a1)); T *= exp(-sig * s); } else T *= exp(-sig * (e1 - e0)); t = (tt); } while (0)
                if (hit && !inner_layers) {
                    ABSORB_TO(b0);
                    C += T; T = 0;
                } else if (hit) {       /* a sphere is one surface to md_render: its front, a layer */
                    ABSORB_TO(b0); C += T * ta; T *= 1 - ta;
                    ABSORB_TO(t1);
                } else {
                    ABSORB_TO(t1);
                }
                #undef ABSORB_TO
                judged++;
                const double e[4] = { C, C, C, 1 - T };
                for (int k = 0; k < 4; ++k) {
                    const double err = fabs(half_bits_to_float(px[4 * p + k]) - e[k]);
                    worst = fmax(worst, err);
                    if (err > 6e-3) { bad++; if (bad < 6) printf("    density%s%s pixel (%u, %u) channel %d: %g against %g\n", ortho ? ", orthographic" : "", inner_layers ? ", layers" : "", x, y, k, half_bits_to_float(px[4 * p + k]), e[k]); }
                }
            }
            printf("    density%s, %s: %d pixels; %.4f off at most\n", ortho ? ", orthographic" : "", inner_layers ? "a transparent sphere's layers inside" : "an opaque sphere inside", judged, worst);
            EXPECT_EQ(0, bad);
            EXPECT_GT(judged, 1000);
        }
    }
    md_gpu_texture_destroy(vc); md_gpu_texture_destroy(vdist); md_gpu_texture_destroy(tft); md_gpu_texture_destroy(out);
    md_compose_destroy(cx);
    free(px); free(fc); free(fd);
    f_close(&fx);
}

/* An isosurface's colour and normal from the volumes beside its field: the
   distance field's sphere (radius 1.3) under BASIC shading without light, so
   a colour shows as it is: times a constant colour volume; through a colour
   map of a colour field (the world's x), with an isoline where it is 0; and,
   lit, with a normals volume of one direction (the view's +z) for every
   point */
UTEST(render, compose_volume_coloring) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 96, H = 72, P = W * H, N = 48, NC = 16;
    const double ctr[3] = { 0, 0, -12 }, R = 1.3;
    cam_t c = { .tan_half = 0.25f, .aspect = (float)W / H, .n = 0.1f };
    md_render_camera_t cam = make_camera(&c);
    md_compose_ctx_t cx = md_compose_create(fx.dev);
    ASSERT_TRUE(cx != NULL);
    md_gpu_texture_t out = f_tex(&fx, MD_GPU_FORMAT_RGBA16_FLOAT, MD_GPU_TEX_STORAGE, W, H, "compose");
    uint16_t* px = malloc(P * 8);
    float* fd = malloc(N * N * N * 4);
    float* fx_ = malloc(NC * NC * NC * 4);
    uint32_t* cv = malloc(NC * NC * NC * 4);
    uint32_t* nv = malloc(NC * NC * NC * 4);
    for (uint32_t z = 0; z < N; ++z) for (uint32_t y = 0; y < N; ++y) for (uint32_t x = 0; x < N; ++x) {
        const double u[3] = { (x + 0.5) / N - 0.5, (y + 0.5) / N - 0.5, (z + 0.5) / N - 0.5 };
        fd[((size_t)z * N + y) * N + x] = (float)(4.0 * sqrt(dot3(u, u)));
    }
    for (uint32_t z = 0; z < NC; ++z) for (uint32_t y = 0; y < NC; ++y) for (uint32_t x = 0; x < NC; ++x) {
        const size_t i = ((size_t)z * NC + y) * NC + x;
        fx_[i] = (float)(4.0 * (x + 0.5) / NC - 2.0);          /* the world's x */
        cv[i]  = 128u | (255u << 8) | (64u << 16) | (255u << 24);
        nv[i]  = 128u | (128u << 8) | (255u << 16) | (255u << 24);   /* +z */
    }
    md_gpu_texture_t vdist = f_tex3d(&fx, N, fd), vx = f_tex3d(&fx, NC, fx_);
    md_gpu_texture_desc_t td = {0};
    td.type = MD_GPU_TEX_3D; td.usage = MD_GPU_TEX_SAMPLED; td.width = td.height = td.depth_or_layers = NC;
    td.format = MD_GPU_FORMAT_RGBA8_UNORM; td.label = "colours";
    md_gpu_texture_t vcol = md_gpu_texture_create(fx.gfx, &td);
    td.label = "normals";
    md_gpu_texture_t vnrm = md_gpu_texture_create(fx.gfx, &td);
    md_gpu_tex_region_t rg = {0};
    ASSERT_TRUE(vdist && vx && vcol && vnrm);
    ASSERT_TRUE(md_gpu_upload_texture(fx.gfx, vcol, &rg, cv, NC * NC * NC * 4) && md_gpu_upload_texture(fx.gfx, vnrm, &rg, nv, NC * NC * NC * 4));
    uint32_t ramp[256];
    for (int i = 0; i < 256; ++i) ramp[i] = (uint32_t)(255 - i) | ((uint32_t)i << 16) | (255u << 24);   /* red to blue */
    md_gpu_texture_t cmap = f_tex(&fx, MD_GPU_FORMAT_RGBA8_UNORM, MD_GPU_TEX_SAMPLED, 256, 1, "colour map");
    ASSERT_TRUE(f_upload_tex(&fx, cmap, ramp, sizeof(ramp)));
    md_render_frame_t f = frame_of(W, H, cam, NULL, 0);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));

    for (int run = 0; run < 3; ++run) {
        md_compose_volume_t v = {0};
        v.field = vdist;
        vol_model(v.model);
        v.clip_max[0] = v.clip_max[1] = v.clip_max[2] = 1.0f;
        v.iso_count = 1; v.iso_values[0] = (float)R;
        for (int k = 0; k < 4; ++k) v.iso_colors[0][k] = 1.0f;
        md_compose_desc_t d = {0};
        d.camera = cam; d.target = out; d.volumes = &v; d.volume_count = 1;
        unlit_shading(&d.shading);
        if (run == 0) v.color_volume = vcol;
        if (run == 1) {
            v.color_field = vx; v.color_map = cmap; v.color_range[0] = (float)-R; v.color_range[1] = (float)R;
            v.isoline_count = 1; v.isoline_values[0] = 0.0f; v.isoline_width = 1.5f;
        }
        md_compose_shading_t* sh = &d.shading;
        if (run == 2) {
            v.normals = vnrm;
            sh->mode = MD_COMPOSE_SHADING_PBR;
            sh->light_dir[0] = 0.4f; sh->light_dir[1] = 0.6f; sh->light_dir[2] = 0.7f; sh->light_strength = 1.0f;
            for (int i = 0; i < 3; ++i) { sh->dir_radiance[i] = 2.7f; sh->sky_radiance[i] = 0.59f; sh->ground_radiance[i] = 0.49f; }
            sh->roughness = 0.4f; sh->f0 = 0.04f;
        }
        ASSERT_TRUE(md_compose(fx.gfx, cx, fx.view, &d));
        ASSERT_TRUE(f_read(&fx, out, px, P * 8));
        int bad = 0, judged = 0, line = 0;
        double worst = 0;
        for (uint32_t y = 0; y < H; ++y) for (uint32_t x = 0; x < W; ++x) {
            const size_t p = (size_t)y * W + x;
            ray_t r = pixel_ray(&c, W, H, x + 0.5, y + 0.5);
            double t0, t1, u0, u1;
            if (!ray_sphere(r, ctr, R - 0.12, &t0, &t1) || !ray_sphere(r, ctr, R, &t0, &t1) || !ray_sphere(r, ctr, R + 0.12, &u0, &u1)) continue;
            const double hx = r.o[0] + r.d[0] * t0;
            const double pxw = 2.0 * c.tan_half * -(r.o[2] + r.d[2] * t0) / H;     /* a pixel's size there */
            double e[3];
            if (run == 0) { e[0] = 128 / 255.0; e[1] = 1.0; e[2] = 64 / 255.0; }
            double tol = 0.012;
            if (run == 1) {
                const double t = (hx + R) / (2 * R);
                /* the map's texels at (i + 0.5) / 256, clamped at its ends */
                const double u = fmin(fmax(t * 256 - 0.5, 0.0), 255.0) / 255.0;
                e[0] = 1 - u; e[1] = 0; e[2] = u;
                /* the line: |x| / (a pixel's x there) pixels from it, near the middle of the view where the
                   surface faces it; dark core over a light rim */
                const double dpx = fabs(hx) / pxw;
                if (dpx < 4.0) {
                    const double w = 1.5, core = fmin(fmax(0.5 * w + 0.5 - dpx, 0.0), 1.0), rim = fmin(fmax(0.5 * (w + 2.0) + 0.5 - dpx, 0.0), 1.0);
                    for (int k = 0; k < 3; ++k) { e[k] += (0.95 - e[k]) * rim; e[k] += (0.02 - e[k]) * core; }
                    tol = 0.06;
                    line += core > 0.5;
                }
            } else if (run == 2) {
                c_lit_t l = c_shade(sh, cf(1, 1, 1), cf((float)-r.d[0], (float)-r.d[1], (float)-r.d[2]), cf(0, 0.0039f, 1), 1.0f);
                e[0] = l.d.x + l.s.x; e[1] = l.d.y + l.s.y; e[2] = l.d.z + l.s.z;
            }
            judged++;
            for (int k = 0; k < 3; ++k) {
                const double err = fabs(half_bits_to_float(px[4 * p + k]) - e[k]);
                worst = fmax(worst, err);
                if (err > tol) { bad++; if (bad < 6) printf("    %s pixel (%u, %u) channel %d: %g against %g\n", run == 0 ? "colour volume" : run == 1 ? "colour field" : "normals", x, y, k, half_bits_to_float(px[4 * p + k]), e[k]); }
            }
        }
        printf("    %s: %d pixels; %.4f off at most%s\n", run == 0 ? "a colour volume" : run == 1 ? "a colour field through a map" : "a normals volume", judged, worst,
               run == 1 ? (line ? ", with the isoline" : ", NO isoline") : "");
        EXPECT_EQ(0, bad);
        EXPECT_GT(judged, 300);
        if (run == 1) EXPECT_GT(line, 10);
    }
    md_gpu_texture_destroy(vdist); md_gpu_texture_destroy(vx); md_gpu_texture_destroy(vcol); md_gpu_texture_destroy(vnrm);
    md_gpu_texture_destroy(cmap); md_gpu_texture_destroy(out);
    md_compose_destroy(cx);
    free(px); free(fd); free(fx_); free(cv); free(nv);
    f_close(&fx);
}

/* Slices: a rectangle x in [-1.5, 1.5], y in [-1, 1] at z -12 facing the
   view (orthographic, a pixel 1/16 of a unit), showing the world's x through
   a red-to-blue map at opacity 0.6: inside it, 0.6 of the map's colour and
   coverage 0.6; over an opaque white sphere behind it, that over the sphere;
   behind an opaque sphere in front of it, the sphere alone. Its edge: within
   its width of the border, white over a dark halo; outside it, nothing. */
UTEST(render, compose_slices) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 128, H = 96, P = W * H, NC = 16;
    cam_t c = { .ortho = true, .half_h = 3.0f, .aspect = (float)W / H, .n = 0.0f, .far_ = 100.0f };
    md_render_camera_t cam = make_camera(&c);
    const double ppu = H / (2.0 * c.half_h);       /* pixels per unit: 16 */
    md_compose_ctx_t cx = md_compose_create(fx.dev);
    ASSERT_TRUE(cx != NULL);
    md_gpu_texture_t out = f_tex(&fx, MD_GPU_FORMAT_RGBA16_FLOAT, MD_GPU_TEX_STORAGE, W, H, "compose");
    uint16_t* px = malloc(P * 8);
    float* fx_ = malloc(NC * NC * NC * 4);
    for (uint32_t z = 0; z < NC; ++z) for (uint32_t y = 0; y < NC; ++y) for (uint32_t x = 0; x < NC; ++x)
        fx_[((size_t)z * NC + y) * NC + x] = (float)(4.0 * (x + 0.5) / NC - 2.0);          /* the world's x over the box */
    md_gpu_texture_t vx = f_tex3d(&fx, NC, fx_);
    uint32_t ramp[256];
    for (int i = 0; i < 256; ++i) ramp[i] = (uint32_t)(255 - i) | ((uint32_t)i << 16) | (255u << 24);
    md_gpu_texture_t cmap = f_tex(&fx, MD_GPU_FORMAT_RGBA8_UNORM, MD_GPU_TEX_SAMPLED, 256, 1, "colour map");
    ASSERT_TRUE(vx && f_upload_tex(&fx, cmap, ramp, sizeof(ramp)));

    /* an opaque white sphere behind the slice's left part, one in front of its right part */
    const float sp[6] = { -0.8f, 0.0f, -16.0f,   0.9f, 0.3f, -9.0f };
    const double sc[2][3] = { { -0.8, 0, -16 }, { 0.9, 0.3, -9 } }, srad[2] = { 0.5, 0.4 };
    float rad[2] = { 0.5f, 0.4f };
    md_render_batch_t b = {0};
    b.key = 1; b.kind = MD_RENDER_SPHERES; b.count = 2; b.points.count = 2; b.points.position = f_buf(&fx, sp, sizeof(sp));
    b.points.radius = f_buf(&fx, rad, sizeof(rad)); b.style.spheres.radius_scale = 1.0f;
    b.color.mode = MD_RENDER_COLOR_UNIFORM; b.color.uniform = 0xFFFFFFFFu;
    md_render_frame_t f = frame_of(W, H, cam, &b, 1);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));

    for (int with_edge = 0; with_edge < 2; ++with_edge) {
        md_compose_slice_t sl = {0};
        memset(sl.rect, 0, sizeof(sl.rect));
        sl.rect[0] = 3.0f; sl.rect[5] = 2.0f; sl.rect[10] = 1.0f; sl.rect[12] = -1.5f; sl.rect[13] = -1.0f; sl.rect[14] = -12.0f; sl.rect[15] = 1.0f;
        /* world to the box's texture space: (x + 2) / 4, (y + 2) / 4, (z + 14) / 4 */
        float* F = sl.field_from_world;
        memset(F, 0, 16 * sizeof(float));
        F[0] = F[5] = F[10] = 0.25f; F[12] = 0.5f; F[13] = 0.5f; F[14] = 3.5f; F[15] = 1.0f;
        sl.field = vx;
        sl.transfer = cmap; sl.transfer_range[0] = -1.5f; sl.transfer_range[1] = 1.5f;
        sl.opacity = 0.6f;
        if (with_edge) { sl.edge[0] = sl.edge[1] = sl.edge[2] = 1.0f; sl.edge[3] = 1.0f; sl.edge_px = 1.0f; }
        md_compose_desc_t d = {0};
        d.camera = cam; d.target = out; d.slices = &sl; d.slice_count = 1;
        unlit_shading(&d.shading);
        ASSERT_TRUE(md_compose(fx.gfx, cx, fx.view, &d));
        ASSERT_TRUE(f_read(&fx, out, px, P * 8));
        int bad = 0, inside = 0, over = 0, under = 0, edge = 0;
        double worst = 0;
        for (uint32_t y = 0; y < H; ++y) for (uint32_t x = 0; x < W; ++x) {
            const size_t p = (size_t)y * W + x;
            ray_t r = pixel_ray(&c, W, H, x + 0.5, y + 0.5);
            const double wx = r.o[0], wy = r.o[1];
            /* pixels in from the rectangle's border (out: negative) */
            const double in_x = fmin(wx + 1.5, 1.5 - wx) * ppu, in_y = fmin(wy + 1.0, 1.0 - wy) * ppu;
            const double dd = (in_x >= 0 && in_y >= 0) ? fmin(in_x, in_y) : -sqrt(fmax(-in_x, 0) * fmax(-in_x, 0) + fmax(-in_y, 0) * fmax(-in_y, 0));
            double t0, t1;
            int hit = -1;
            for (int k = 0; k < 2; ++k) {
                double u0, u1;
                if (ray_sphere(r, sc[k], srad[k] - 0.1, &u0, &u1) != ray_sphere(r, sc[k], srad[k] + 0.1, &u0, &u1)) hit = -2;
                else if (hit == -1 && ray_sphere(r, sc[k], srad[k], &t0, &t1)) hit = k;
            }
            if (hit == -2) continue;
            const double got[4] = { half_bits_to_float(px[4 * p]), half_bits_to_float(px[4 * p + 1]), half_bits_to_float(px[4 * p + 2]), half_bits_to_float(px[4 * p + 3]) };
            if (with_edge) {
                if (hit >= 0) continue;
                if (dd < -3.0) { for (int k = 0; k < 4; ++k) bad += got[k] != 0.0; continue; }     /* past the edge and its halo */
                if (fabs(dd) <= 0.51) {         /* on the line, half a pixel from the border: its core, white, all of it */
                    edge++;
                    for (int k = 0; k < 4; ++k) if (fabs(got[k] - 1.0) > 2e-3) { bad++; if (bad < 6) printf("    edge pixel (%u, %u) channel %d: %g\n", x, y, k, got[k]); }
                }
                continue;
            }
            if (fabs(dd) < 0.6) continue;
            double e[4] = { 0, 0, 0, 0 };
            const double al = dd > 0 ? 0.6 : 0.0;
            if (al > 0 && hit != 1) {
                const double t = (wx + 1.5) / 3.0, u = fmin(fmax(t * 256 - 0.5, 0.0), 255.0) / 255.0;
                e[0] = al * c_s2l((float)(1 - u)); e[2] = al * c_s2l((float)u); e[3] = al;
                inside++;
            }
            if (hit == 1) { for (int k = 0; k < 4; ++k) e[k] = 1.0; under++; }                  /* in front of it */
            else if (hit == 0) { for (int k = 0; k < 3; ++k) e[k] += (1 - al); e[3] = 1.0; over += al > 0; }
            for (int k = 0; k < 4; ++k) {
                const double err = fabs(got[k] - e[k]);
                worst = fmax(worst, err);
                if (err > 6e-3) { bad++; if (bad < 6) printf("    slice pixel (%u, %u) channel %d: %g against %g\n", x, y, k, got[k], e[k]); }
            }
        }
        if (with_edge) {
            printf("    its edge: %d pixels on it white, nothing past it\n", edge);
            EXPECT_GT(edge, 100);
        } else {
            printf("    a slice: %d pixels on it (%d over a sphere behind it), %d behind a sphere in front of it; %.4f off at most\n", inside, over, under, worst);
            EXPECT_GT(inside, 1000);
            EXPECT_GT(over, 50);
            EXPECT_GT(under, 50);
        }
        EXPECT_EQ(0, bad);
    }
    md_gpu_texture_destroy(vx); md_gpu_texture_destroy(cmap); md_gpu_texture_destroy(out);
    md_compose_destroy(cx);
    free(px); free(fx_);
    f_close(&fx);
}

/* ---- Post ----------------------------------------------------------------------------- */

/* md_post.slang's display transform, in C: the reference */
static cf3 c_neutral(cf3 c) {
    const float sc = 0.8f - 0.04f, des = 0.15f;
    float x = fminf(c.x, fminf(c.y, c.z));
    float off = x < 0.08f ? x - 6.25f * x * x : 0.04f;
    c = cf(c.x - off, c.y - off, c.z - off);
    float peak = fmaxf(c.x, fmaxf(c.y, c.z));
    if (peak < sc) return c;
    const float d = 1.0f - sc;
    float np = 1.0f - d * d / (peak + d - sc);
    c = cmul(c, np / peak);
    float g = 1.0f - 1.0f / (des * (peak - np) + 1.0f);
    return cmix3(c, cf(np, np, np), g);
}
static cf3 c_aces(cf3 c) {
    cf3 v = cf(0.59719f * c.x + 0.35458f * c.y + 0.04823f * c.z, 0.07600f * c.x + 0.90834f * c.y + 0.01566f * c.z, 0.02840f * c.x + 0.13383f * c.y + 0.83777f * c.z);
    float r[3], vv[3] = { v.x, v.y, v.z };
    for (int i = 0; i < 3; ++i) r[i] = (vv[i] * (vv[i] + 0.0245786f) - 0.000090537f) / (vv[i] * (0.983729f * vv[i] + 0.4329510f) + 0.238081f);
    return cf(1.60475f * r[0] - 0.53108f * r[1] - 0.07367f * r[2], -0.10208f * r[0] + 1.10813f * r[1] - 0.00605f * r[2], -0.00327f * r[0] - 0.07276f * r[1] + 1.07602f * r[2]);
}
static float c_hable1(float x) { return ((x * (0.15f * x + 0.05f) + 0.004f) / (x * (0.15f * x + 0.5f) + 0.06f)) - 0.02f / 0.3f; }
static cf3 c_curve(int mode, cf3 x) {
    if (mode == 1) { float w = c_neutral(cf(8, 8, 8)).x; return cmul(c_neutral(x), 1.0f / w); }
    if (mode == 2) { float w = c_aces(cf(8, 8, 8)).x;    return cmul(c_aces(x), 1.0f / w); }
    if (mode == 3) { float w = c_hable1(11.2f); return cf(c_hable1(2 * x.x) / w, c_hable1(2 * x.y) / w, c_hable1(2 * x.z) / w); }
    return x;
}
static float c_l2s(float c) { c = csat(c); return c < 0.0031308f ? c * 12.92f : 1.055f * powf(c, 1.0f / 2.4f) - 0.055f; }

static md_gpu_texture_t f_tex(fix_t* fx, md_gpu_format_t format, md_gpu_tex_usage_t usage, uint32_t w, uint32_t h, const char* label) {
    md_gpu_texture_desc_t td = {0};
    td.type = MD_GPU_TEX_2D; td.usage = usage; td.width = w; td.height = h; td.format = format; td.label = label;
    return md_gpu_texture_create(fx->gfx, &td);
}

static bool f_upload_tex(fix_t* fx, md_gpu_texture_t t, const void* src, size_t size) {
    md_gpu_tex_region_t r = {0};
    r.extent[2] = 1;
    return md_gpu_upload_texture(fx->gfx, t, &r, src, size);
}

/* The display transform of every curve over premultiplied radiance of
   coverage 0, a quarter, 0.6 and 1, against the reference */
UTEST(render, post_display_curves) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 32, H = 8, P = W * H;
    uint16_t* in = malloc(P * 8);
    uint16_t* px = malloc(P * 8);
    for (uint32_t i = 0; i < P; ++i) {
        const float a = (float[]){ 0.0f, 0.25f, 0.6f, 1.0f }[i % 4];
        const float L[3] = { 12.0f * fmodf(i * 0.6180339f, 1.0f), 3.0f * fmodf(i * 0.4142135f, 1.0f), 0.5f * fmodf(i * 0.7320508f, 1.0f) };
        for (int k = 0; k < 3; ++k) in[4 * i + k] = float_to_half(L[k] * a);
        in[4 * i + 3] = float_to_half(a);
    }
    md_gpu_texture_t scene = f_tex(&fx, MD_GPU_FORMAT_RGBA16_FLOAT, MD_GPU_TEX_SAMPLED, W, H, "scene");
    md_gpu_texture_t out   = f_tex(&fx, MD_GPU_FORMAT_RGBA16_FLOAT, MD_GPU_TEX_STORAGE, W, H, "display");
    ASSERT_TRUE(f_upload_tex(&fx, scene, in, P * 8));
    md_post_ctx_t pc = md_post_create(fx.dev);
    ASSERT_TRUE(pc != NULL);
    for (int curve = 0; curve < 4; ++curve) {
        md_post_display_t d = { .scene = scene, .target = out, .curve = (md_post_curve_t)curve, .exposure = 1.3f, .gamma = 2.0f,
                                .background = { 0.9f, 0.5f, 0.1f } };
        ASSERT_TRUE(md_post_display(fx.gfx, pc, &d));
        ASSERT_TRUE(f_read(&fx, out, px, P * 8));
        /* without a curve: exposure 1, gamma 2.2 */
        const float e = curve ? 1.3f : 1.0f, g = curve ? 2.0f : 2.2f, kg = g / 2.2f;
        const cf3 bg = cf(c_s2l(powf(0.9f, kg)), c_s2l(powf(0.5f, kg)), c_s2l(powf(0.1f, kg)));
        double worst = 0;
        int bad = 0;
        for (uint32_t i = 0; i < P; ++i) {
            const float a = half_bits_to_float(in[4 * i + 3]);
            cf3 L = cf(half_bits_to_float(in[4 * i]), half_bits_to_float(in[4 * i + 1]), half_bits_to_float(in[4 * i + 2]));
            cf3 sc = a > 1e-4f ? c_curve(curve, cmul(L, e / a)) : cf(0, 0, 0);
            cf3 col = cadd(cmul(cf(csat(sc.x), csat(sc.y), csat(sc.z)), a), cmul(bg, 1 - a));
            const float exp[3] = { powf(c_l2s(col.x), 2.2f / g), powf(c_l2s(col.y), 2.2f / g), powf(c_l2s(col.z), 2.2f / g) };
            for (int k = 0; k < 3; ++k) {
                const double err = fabs(half_bits_to_float(px[4 * i + k]) - exp[k]);
                worst = fmax(worst, err);
                if (err > 2e-3) { bad++; if (bad < 6) printf("    curve %d pixel %u channel %d: %g against %g\n", curve, i, k, half_bits_to_float(px[4 * i + k]), exp[k]); }
            }
            if (half_bits_to_float(px[4 * i + 3]) != 1.0f) bad++;
        }
        printf("    curve %d: %.5f off at most\n", curve, worst);
        EXPECT_EQ(0, bad);
    }
    md_post_destroy(pc);
    md_gpu_texture_destroy(scene);
    md_gpu_texture_destroy(out);
    free(in); free(px);
    f_close(&fx);
}

/* A view's frame through md_compose and md_post_display */
typedef struct {
    md_compose_ctx_t cx;
    md_post_ctx_t    pc;
    md_gpu_texture_t scene, frame;
} f_post_t;

static bool f_frame(fix_t* fx, f_post_t* fp, md_render_frame_t* f, md_render_status_t* st) {
    if (!md_render_view_draw(fx->gfx, fx->view, f, st)) return false;
    md_compose_desc_t d = {0};
    d.camera = f->camera;
    d.target = fp->scene;
    d.shading.mode = MD_COMPOSE_SHADING_PBR;
    d.shading.light_dir[0] = 0.4f; d.shading.light_dir[1] = 0.6f; d.shading.light_dir[2] = 0.7f;
    d.shading.light_strength = 1.0f;
    for (int i = 0; i < 3; ++i) { d.shading.dir_radiance[i] = 2.7f; d.shading.sky_radiance[i] = 0.59f; d.shading.ground_radiance[i] = 0.49f; }
    d.shading.roughness = 0.4f; d.shading.f0 = 0.04f;
    if (!md_compose(fx->gfx, fp->cx, fx->view, &d)) return false;
    md_post_display_t pd = { .scene = fp->scene, .target = fp->frame, .curve = MD_POST_CURVE_NEUTRAL, .exposure = 1.0f, .gamma = 2.2f,
                             .background = { 1, 1, 1 } };
    return md_post_display(fx->gfx, fp->pc, &pd);
}

/* The temporal resolve. A frame jittered by j, resolved alone, is the
   unjittered frame (bilinearly, where it is smooth), and far less so with
   the jitter's sign turned; a still view resolved over 16 jittered frames
   settles on the unjittered image. */
UTEST(render, post_resolve_unjitters_and_settles) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 96, H = 72, P = W * H;
    cam_t c = { .tan_half = 0.4f, .aspect = (float)W / H, .n = 0.1f };
    const float ctr[3] = { 0.0f, 0.0f, -12.0f };
    md_render_batch_t b = {0};
    b.key = 1; b.kind = MD_RENDER_SPHERES; b.count = 1; b.points.count = 1; b.points.position = f_buf(&fx, ctr, sizeof(ctr));
    b.style.spheres.radius_add = 3.0f; b.color.mode = MD_RENDER_COLOR_UNIFORM; b.color.uniform = 0xFF3080D0u;
    f_post_t fp = { md_compose_create(fx.dev), md_post_create(fx.dev),
                    f_tex(&fx, MD_GPU_FORMAT_RGBA16_FLOAT, MD_GPU_TEX_STORAGE | MD_GPU_TEX_SAMPLED, W, H, "scene"),
                    f_tex(&fx, MD_GPU_FORMAT_RGBA16_FLOAT, MD_GPU_TEX_STORAGE | MD_GPU_TEX_SAMPLED, W, H, "frame") };
    ASSERT_TRUE(fp.cx && fp.pc);
    md_gpu_texture_t hist[2] = { f_tex(&fx, MD_GPU_FORMAT_RGBA16_FLOAT, MD_GPU_TEX_STORAGE | MD_GPU_TEX_SAMPLED, W, H, "history 0"),
                                 f_tex(&fx, MD_GPU_FORMAT_RGBA16_FLOAT, MD_GPU_TEX_STORAGE | MD_GPU_TEX_SAMPLED, W, H, "history 1") };
    uint16_t* ref = malloc(P * 8);
    uint16_t* got = malloc(P * 8);
    uint16_t* last = malloc(P * 8);

    /* the unjittered frame */
    md_render_frame_t f = frame_of(W, H, make_camera(&c), &b, 1);
    f.jitter = MD_RENDER_JITTER_NONE;
    ASSERT_TRUE(f_frame(&fx, &fp, &f, NULL));
    ASSERT_TRUE(f_read(&fx, fp.frame, ref, P * 8));
    /* the sphere's inside, away from its edge: r < 0.7 of its radius on the screen */
    const double rs = 3.0 / (12.0 * c.tan_half) * 0.5 * H;
    #define INSIDE(x, y) (((x) + 0.5 - 0.5 * W) * ((x) + 0.5 - 0.5 * W) + ((y) + 0.5 - 0.5 * H) * ((y) + 0.5 - 0.5 * H) < 0.49 * rs * rs)
    double err[2] = { 0, 0 }, worst = 0;
    int n = 0;
    const float j[2] = { 0.3f, -0.2f };
    f.jitter = MD_RENDER_JITTER_GIVEN;
    f.jitter_offset[0] = j[0]; f.jitter_offset[1] = j[1];
    ASSERT_TRUE(f_frame(&fx, &fp, &f, NULL));
    for (int sign = 0; sign < 2; ++sign) {
        md_post_resolve_t r = { .view = fx.view, .camera = f.camera, .jitter = { sign ? -j[0] : j[0], sign ? -j[1] : j[1] },
                                .frame = fp.frame, .history = hist[0], .feedback_min = 0.8f, .feedback_max = 0.95f };
        ASSERT_TRUE(md_post_resolve(fx.gfx, fp.pc, &r));
        ASSERT_TRUE(f_read(&fx, hist[0], got, P * 8));
        n = 0;
        for (uint32_t y = 0; y < H; ++y) for (uint32_t x = 0; x < W; ++x) {
            if (!INSIDE(x, y)) continue;
            const size_t p = (size_t)y * W + x;
            for (int k = 0; k < 3; ++k) {
                const double e = fabs(half_bits_to_float(got[4 * p + k]) - half_bits_to_float(ref[4 * p + k]));
                err[sign] += e;
                if (!sign) worst = fmax(worst, e);
            }
            n++;
        }
        err[sign] /= 3.0 * n;
    }
    printf("    a frame jittered by (%.1f, %.1f) resolved alone, against the unjittered over %d pixels: %.5f mean (%.5f at most); the jitter's sign turned: %.5f\n",
           j[0], j[1], n, err[0], worst, err[1]);
    EXPECT_LT(err[0], 0.2 * err[1]);
    EXPECT_LT(err[0], 0.002);
    EXPECT_LT(worst, 0.03);     /* bilinear across the highlight's peak */

    /* still: 16 frames of the Halton sequence, each resolved into the other history */
    int cur = 0, nan = 0;
    double settle = 0, still = 0;
    for (int i = 0; i < 16; ++i) {
        float hx = 0, hy = 0;
        for (int base = 2, k = 0; k < 2; base = 3, ++k) {
            float r = 0, fr = 1.0f / base;
            for (int m = i + 1; m > 0; m /= base, fr /= base) r += fr * (m % base);
            if (k == 0) hx = r - 0.5f; else hy = r - 0.5f;
        }
        f.jitter_offset[0] = hx; f.jitter_offset[1] = hy;
        ASSERT_TRUE(f_frame(&fx, &fp, &f, NULL));
        md_post_resolve_t r = { .view = fx.view, .camera = f.camera, .jitter = { hx, hy }, .frame = fp.frame,
                                .history_prev = i ? hist[cur ^ 1] : (md_gpu_texture_t)0, .history = hist[cur], .feedback_min = 0.8f, .feedback_max = 0.95f };
        ASSERT_TRUE(md_post_resolve(fx.gfx, fp.pc, &r));
        if (i == 14) ASSERT_TRUE(f_read(&fx, hist[cur], last, P * 8));
        cur ^= 1;
    }
    ASSERT_TRUE(f_read(&fx, hist[cur ^ 1], got, P * 8));
    for (uint32_t p = 0; p < P; ++p) for (int k = 0; k < 4; ++k) nan += !isfinite(half_bits_to_float(got[4 * p + k]));
    for (uint32_t y = 0; y < H; ++y) for (uint32_t x = 0; x < W; ++x) {
        if (!INSIDE(x, y)) continue;
        const size_t p = (size_t)y * W + x;
        for (int k = 0; k < 3; ++k) {
            settle = fmax(settle, fabs(half_bits_to_float(got[4 * p + k]) - half_bits_to_float(ref[4 * p + k])));
            still  = fmax(still,  fabs(half_bits_to_float(got[4 * p + k]) - half_bits_to_float(last[4 * p + k])));
        }
    }
    #undef INSIDE
    printf("    still, 16 jittered frames: %.5f off the unjittered at most, %.5f from the frame before\n", settle, still);
    EXPECT_EQ(0, nan);
    EXPECT_LT(settle, 0.03);
    EXPECT_LT(still, 0.005);
    md_gpu_texture_destroy(hist[0]); md_gpu_texture_destroy(hist[1]);
    md_gpu_texture_destroy(fp.scene); md_gpu_texture_destroy(fp.frame);
    md_post_destroy(fp.pc); md_compose_destroy(fp.cx);
    free(ref); free(got); free(last);
    f_close(&fx);
}

/* The marks over an image: a selected sphere in front of a highlighted one
   of tag 3 (its group's colour): the selection's tint where it is, the
   highlight's in its group's colour where it shows, its fill where it is
   behind the selected sphere, the image (dithered) elsewhere */
UTEST(render, post_present_marks) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 96, H = 64, P = W * H;
    cam_t c = { .ortho = true, .half_h = 4.0f, .aspect = (float)W / H, .n = 0.0f, .far_ = 100.0f };
    float pos[6] = { -1.0f, 0.0f, -10.0f,   1.0f, 0.0f, -20.0f };
    uint8_t marks[2] = { MD_RENDER_MARK_STATE0, MD_RENDER_MARK_STATE1 | MD_RENDER_MARK_TAG(3) };
    md_render_batch_t b = {0};
    b.key = 1; b.kind = MD_RENDER_SPHERES; b.count = 2; b.points.count = 2;
    b.points.position = f_buf(&fx, pos, sizeof(pos)); b.points.marks = f_buf(&fx, marks, sizeof(marks));
    b.style.spheres.radius_add = 1.6f;
    md_render_frame_t f = frame_of(W, H, make_camera(&c), &b, 1);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));

    uint16_t* img = malloc(P * 8);
    for (uint32_t i = 0; i < P * 4; ++i) img[i] = float_to_half((i & 3) == 3 ? 1.0f : 0.5f);
    md_gpu_texture_t image  = f_tex(&fx, MD_GPU_FORMAT_RGBA16_FLOAT, MD_GPU_TEX_SAMPLED, W, H, "image");
    md_gpu_texture_t target = f_tex(&fx, MD_GPU_FORMAT_RGBA8_UNORM, MD_GPU_TEX_STORAGE, W, H, "shown");
    ASSERT_TRUE(f_upload_tex(&fx, image, img, P * 8));
    md_post_marks_t mk = {0};
    mk.visible[0][0] = 1.0f; mk.visible[0][3] = 0.5f;                       /* the selection: red, half */
    mk.visible[1][2] = 1.0f; mk.visible[1][3] = 0.5f;                       /* the highlight: blue, half (but its group's) */
    mk.hidden[1][1]  = 1.0f; mk.hidden[1][3]  = 1.0f;                       /* hidden: green, all of it (but its group's) */
    mk.groups[2][0] = 1.0f; mk.groups[2][1] = 1.0f;                         /* tag 3: yellow */
    mk.highlight = 1.0f;
    mk.outline_px = 1;
    md_post_ctx_t pc = md_post_create(fx.dev);
    ASSERT_TRUE(pc != NULL);
    md_post_present_t pr = { .image = image, .marks_view = fx.view, .marks = &mk, .target = target, .frame = 7 };
    ASSERT_TRUE(md_post_present(fx.gfx, pc, &pr));
    uint32_t* out = malloc(P * 4);
    ASSERT_TRUE(f_read(&fx, target, out, P * 4));
    /* at a world point: the pixel's channel k, 0..1 */
    #define PIX(wx, wy, k) (((out[(uint32_t)((0.5 - (wy) / (2.0 * c.half_h)) * H) * W + (uint32_t)(((wx) / (2.0 * c.half_h * c.aspect)) * W + 0.5 * W)] >> (8 * (k))) & 0xFFu) / 255.0)
    const struct { float x, y; float rgb[3]; const char* what; } at[] = {
        { -1.8f, 0.0f, { 0.75f, 0.25f, 0.25f }, "the selection" },
        {  1.8f, 0.0f, { 0.75f, 0.75f, 0.25f }, "the highlight in view (its group's)" },
        {  0.3f, 0.0f, { 1.00f, 1.00f, 0.00f }, "the highlight behind the selection (its fill, its group's)" },
        {  3.8f, 2.5f, { 0.50f, 0.50f, 0.50f }, "nothing" },
    };
    int bad = 0;
    for (size_t i = 0; i < sizeof(at) / sizeof(at[0]); ++i) {
        for (int k = 0; k < 3; ++k) {
            const double v = PIX(at[i].x, at[i].y, k);
            if (fabs(v - at[i].rgb[k]) > 2.5 / 255.0) { bad++; printf("    %s: channel %d %.3f against %.3f\n", at[i].what, k, v, at[i].rgb[k]); }
        }
    }
    #undef PIX
    EXPECT_EQ(0, bad);
    /* without marks: the image, dithered by at most a step */
    md_post_present_t pn = { .image = image, .target = target, .frame = 8 };
    ASSERT_TRUE(md_post_present(fx.gfx, pc, &pn));
    ASSERT_TRUE(f_read(&fx, target, out, P * 4));
    int off = 0;
    for (uint32_t p = 0; p < P; ++p) for (int k = 0; k < 3; ++k) off += abs((int)((out[p] >> (8 * k)) & 0xFFu) - 128) > 1;
    EXPECT_EQ(0, off);
    md_post_destroy(pc);
    md_gpu_texture_destroy(image); md_gpu_texture_destroy(target);
    free(img); free(out);
    f_close(&fx);
}

/* md_render.h's side of a displaced capsule, in double: pa, pb in view space,
   side (already in view space, divided by the scale) or NULL */
static void ref_capsule_side(const double pa[3], const double pb[3], const double* sv, bool ortho, double out[3]) {
    double z[3] = { pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2] }, l = sqrt(dot3(z, z));
    for (int k = 0; k < 3; ++k) z[k] /= l;
    double w[3] = { 0, 0, 1 };
    if (!ortho) for (int k = 0; k < 3; ++k) w[k] = -0.5 * (pa[k] + pb[k]);
    double u[3], vv[3];
    cross3(z, w, u);
    double lu = sqrt(dot3(u, u));
    for (int k = 0; k < 3; ++k) u[k] /= lu;
    if (!sv) { memcpy(out, u, sizeof(u)); return; }
    cross3(z, u, vv);
    double x = dot3(sv, u), y = dot3(sv, vv), m2 = x * x + y * y;
    double q[2] = { 0.02, 0.0 };
    if (m2 > 0) { q[0] += (x * x - y * y) / sqrt(m2); q[1] += 2 * x * y / sqrt(m2); }
    double ang = 0.5 * atan2(q[1], q[0]);         /* half the angle */
    for (int k = 0; k < 3; ++k) out[k] = cos(ang) * u[k] + sin(ang) * vv[k];
    if (dot3(out, sv) < 0) for (int k = 0; k < 3; ++k) out[k] = -out[k];
}

/* Displaced capsules (md_render.h: offset, side) against plain capsules
   moved by the reference's side: strong sides, weak ones, none, one with a
   part along the axis; offsets of either sign over the same capsules (a
   double bond), and without sides (facing the view); a second instance
   turned and scaled, its sides weighed alike */
UTEST(render, capsules_displaced_follow_their_sides) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 192, H = 144;
    cam_t c = { .tan_half = 0.4f, .aspect = (float)W / H, .n = 0.1f };
    md_render_camera_t cam = make_camera(&c);
    enum { NC = 6 };
    const float pos[2 * NC][3] = {
        { -5.0f,  3.0f, -14.0f }, { -1.0f,  3.2f, -14.0f },
        {  0.5f,  3.0f, -13.0f }, {  4.5f,  1.5f, -16.0f },
        { -5.0f,  0.0f, -12.0f }, { -2.6f,  0.9f, -16.5f },
        {  0.0f,  0.0f, -14.0f }, {  3.6f,  0.4f, -14.0f },
        { -5.5f, -2.5f, -14.0f }, { -1.5f, -3.5f, -15.0f },
        {  1.0f, -2.0f, -14.0f }, {  4.4f, -3.4f, -15.0f },
    };
    const float side[NC][4] = {
        { 0.0f, 1.0f, 0.0f, 0 },          /* strong, across the view: holds */
        { 0.0f, 0.0f, 1.0f, 0 },          /* strong, along the view: holds, sticks one behind the other */
        { 0.0f, 0.004f, 0.01f, 0 },       /* weak: mostly the view's */
        { 0.0f, 0.0f, 0.0f, 0 },          /* none: the view's */
        { 0.6f, 0.4f, 0.7f, 0 },          /* with a part along the axis */
        { 0.05f, 0.05f, 0.05f, 0 },       /* middling */
    };
    const float R = 0.16f, OFF = 1.38f, OFF3 = 2.57f;
    float xf[32] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };
    {
        const double ang = 35 * 3.14159265358979 / 180, sc = 1.3, ctr[3] = { 0, 0, -14 }, mv[3] = { 0.3, -0.5, -7.0 };
        float m[9] = { (float)(sc * cos(ang)), 0, (float)(-sc * sin(ang)), 0, (float)sc, 0, (float)(sc * sin(ang)), 0, (float)(sc * cos(ang)) };
        for (int i = 0; i < 9; ++i) xf[16 + 4 * (i / 3) + i % 3] = m[i];
        for (int i = 0; i < 3; ++i) xf[28 + i] = (float)(ctr[i] + mv[i] - (m[i] * ctr[0] + m[3 + i] * ctr[1] + m[6 + i] * ctr[2]));
    }
    /* a fourth: pulled in by miters and a trim (the partners of delocalized bonds) */
    float side4[NC][4];
    const float miter[NC][2] = { { 0.577f, 0.577f }, { 0.0f, 1.2f }, { 0.4f, 0.0f }, { 0.7f, 0.3f }, { 2.0f, 2.5f }, { 0.0f, 0.0f } };
    for (int k = 0; k < NC; ++k) {
        memcpy(side4[k], side[k], sizeof(side4[k]));
        const uint32_t w = (uint32_t)float_to_half(miter[k][0]) | ((uint32_t)float_to_half(miter[k][1]) << 16);
        memcpy(&side4[k][3], &w, 4);
    }
    const float OFF4 = 3.0f, TRIM = 0.08f;
    md_render_batch_t b[4];
    memset(b, 0, sizeof(b));
    b[0].key = 51; b[0].kind = MD_RENDER_CAPSULES; b[0].count = NC; b[0].points.count = 2 * NC;
    b[0].points.position = f_buf(&fx, pos, sizeof(pos));
    b[0].style.capsules.radius = R;
    b[0].style.capsules.offset = OFF;
    b[0].style.capsules.side   = f_buf(&fx, side, sizeof(side));
    b[0].transform = xf; b[0].transform_count = 2;
    b[1] = b[0]; b[1].key = 52; b[1].style.capsules.offset = -OFF;
    b[2] = b[0]; b[2].key = 53; b[2].style.capsules.offset = OFF3; b[2].style.capsules.side = 0;
    b[3] = b[0]; b[3].key = 54; b[3].style.capsules.offset = -OFF4; b[3].style.capsules.side = f_buf(&fx, side4, sizeof(side4));
    b[3].style.capsules.trim = TRIM;
    md_render_frame_t f = frame_of(W, H, cam, b, 4);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));

    double ca[4 * 2 * NC][3], cb[4 * 2 * NC][3], cr[4 * 2 * NC];
    uint32_t ckey[4 * 2 * NC], clin[4 * 2 * NC];
    int n = 0;
    for (int bi = 0; bi < 4; ++bi) {
        for (int in = 0; in < 2; ++in) {
            const float* m = xf + 16 * in;
            const double sc = sqrt((double)m[0] * m[0] + (double)m[1] * m[1] + (double)m[2] * m[2]);
            for (int k = 0; k < NC; ++k) {
                double p[2][3], sv[3];
                for (int e = 0; e < 2; ++e) for (int i = 0; i < 3; ++i) {
                    const float* q = pos[2 * k + e];
                    p[e][i] = (double)m[i] * q[0] + (double)m[4 + i] * q[1] + (double)m[8 + i] * q[2] + m[12 + i];
                }
                for (int i = 0; i < 3; ++i) sv[i] = ((double)m[i] * side[k][0] + (double)m[4 + i] * side[k][1] + (double)m[8 + i] * side[k][2]) / sc;
                double d[3];
                ref_capsule_side(p[0], p[1], bi != 2 ? sv : NULL, false, d);
                const double o = (bi == 0 ? OFF : bi == 1 ? -OFF : bi == 2 ? OFF3 : -OFF4) * R * sc;
                for (int i = 0; i < 3; ++i) { ca[n][i] = p[0][i] + d[i] * o; cb[n][i] = p[1][i] + d[i] * o; }
                if (bi == 3) {          /* pulled in: by the miters times |o|, and the trim, to a fifth at least */
                    double ax[3] = { cb[n][0] - ca[n][0], cb[n][1] - ca[n][1], cb[n][2] - ca[n][2] }, len = sqrt(dot3(ax, ax));
                    double ia = fabs(o) * miter[k][0] + TRIM * len, ib = fabs(o) * miter[k][1] + TRIM * len;
                    if (ia + ib > 0.8 * len) { const double q = 0.8 * len / (ia + ib); ia *= q; ib *= q; }
                    for (int i = 0; i < 3; ++i) { ca[n][i] += ax[i] / len * ia; cb[n][i] -= ax[i] / len * ib; }
                }
                cr[n] = R * sc; ckey[n] = 51u + (uint32_t)bi; clin[n] = (uint32_t)(in * NC + k);
                n++;
            }
        }
    }
    md_render_targets_t tg = md_render_view_targets(fx.view);
    uint32_t* id = malloc(W * H * 8);
    float*    dp = malloc(W * H * 4);
    ASSERT_TRUE(f_read(&fx, tg.id, id, W * H * 8) && f_read(&fx, tg.depth, dp, W * H * 4));
    int bad = 0, covered = 0, judged = 0;
    double worst = 0;
    for (uint32_t y = 0; y < H; ++y) {
        for (uint32_t x = 0; x < W; ++x) {
            const uint32_t i = y * W + x;
            ray_t r = pixel_ray(&c, W, H, x + 0.5, y + 0.5);
            const double far_[3] = { r.d[0] * 100, r.d[1] * 100, r.d[2] * 100 };
            bool hit_in = false, hit_out = false;
            for (int k = 0; k < n; ++k) {
                const double d = seg_seg_dist(r.o, far_, ca[k], cb[k]);
                const double px = -0.5 * (ca[k][2] + cb[k][2]) * 2.0 * c.tan_half / H;
                hit_in  |= d < cr[k] - 0.6 * px;
                hit_out |= d < cr[k] + 0.6 * px;
            }
            if (id[2 * i + 1] == 0) {
                if (hit_in) { bad++; if (bad < 8) printf("    pixel (%u, %u): empty, a capsule there\n", x, y); }
                continue;
            }
            covered++;
            if (!hit_out) { bad++; if (bad < 8) printf("    pixel (%u, %u): covered, no capsule near\n", x, y); continue; }
            const double zv = -c.n / dp[i], tp = zv / r.d[2];
            const double P[3] = { r.o[0] + r.d[0] * tp, r.o[1] + r.d[1] * tp, r.o[2] + r.d[2] * tp };
            const double pxw = -P[2] * 2.0 * c.tan_half / H;
            int who = -1;
            double best = 1e30, second = 1e30;
            for (int k = 0; k < n; ++k) {
                const double sd = pt_seg_dist(P, ca[k], cb[k]) - cr[k];
                if (sd < best) { second = best; best = sd; who = k; } else if (sd < second) second = sd;
            }
            worst = fmax(worst, fabs(best) / pxw);
            if (fabs(best) > 0.05 * pxw) { bad++; if (bad < 8) printf("    pixel (%u, %u): %.3f px off the capsules\n", x, y, best / pxw); continue; }
            if (second - best > pxw) {
                judged++;
                const uint32_t slot = (id[2 * i + 1] >> 20) - 1u;
                if (id[2 * i] != clin[who] || slot != ckey[who] - 51u) {
                    bad++;
                    if (bad < 8) printf("    pixel (%u, %u): batch %u element %u, expected %u %u\n", x, y, slot, id[2 * i], ckey[who] - 51u, clin[who]);
                }
            }
        }
    }
    printf("    %d pixels covered, %d of one capsule clearly; the farthest %.4f px from the capsules\n", covered, judged, worst);
    EXPECT_EQ(0, bad);
    EXPECT_GT(covered, 2000);
    EXPECT_GT(judged, covered / 2);

    /* the strong side across the view holds: the two sticks of the first, one above the other */
    EXPECT_NEAR(ca[0][1] - pos[0][1], OFF * R, 0.005 * OFF * R);      /* the view's 1/50 turns it a little */
    EXPECT_NEAR(ca[2 * NC][1] - pos[0][1], -OFF * R, 0.005 * OFF * R);
    /* a pick on the middle of the upper stick of the first: that batch and element */
    {
        double q[3];
        for (int i = 0; i < 3; ++i) q[i] = 0.5 * (ca[0][i] + cb[0][i]);
        const uint32_t qx = (uint32_t)((q[0] / -q[2] / (c.tan_half * c.aspect) * 0.5 + 0.5) * W);
        const uint32_t qy = (uint32_t)((0.5 - q[1] / -q[2] / c.tan_half * 0.5) * H);
        md_render_pick_t pk;
        ASSERT_TRUE(f_pick(&fx, qx, qy, &pk));
        EXPECT_EQ(51u, (unsigned)pk.key);
        EXPECT_EQ(0u, pk.element);
    }
    free(id); free(dp);
    f_close(&fx);
}

/* A circle of md_render.h (POLYGONS with an inset), in double: its frame,
   its radii, its corners' angles */
typedef struct {
    double c[3], n[3], t[3], b[3], R, r, ang[8];
    int    nc, dashes;
    uint32_t key, linear;
} circ_ref_t;

static bool circ_ref_build(circ_ref_t* o, const double (*q)[3], int n, double inset, double r, int dashes) {
    memset(o, 0, sizeof(*o));
    o->nc = n; o->r = r; o->dashes = dashes;
    for (int k = 0; k < n; ++k) for (int i = 0; i < 3; ++i) o->c[i] += q[k][i] / n;
    double nv[3] = { 0, 0, 0 };
    for (int k = 0; k < n; ++k) {
        double a[3], b[3], x[3];
        for (int i = 0; i < 3; ++i) { a[i] = q[k][i] - o->c[i]; b[i] = q[(k + 1) % n][i] - o->c[i]; }
        cross3(a, b, x);
        for (int i = 0; i < 3; ++i) nv[i] += x[i];
    }
    double nl = sqrt(dot3(nv, nv));
    for (int i = 0; i < 3; ++i) o->n[i] = nv[i] / nl;
    double p[8][3];
    for (int k = 0; k < n; ++k) {
        double d[3] = { q[k][0] - o->c[0], q[k][1] - o->c[1], q[k][2] - o->c[2] }, h = dot3(d, o->n);
        for (int i = 0; i < 3; ++i) p[k][i] = d[i] - o->n[i] * h;
    }
    double ap = 1e30;
    for (int k = 0; k < n; ++k) {
        double d[3] = { p[(k + 1) % n][0] - p[k][0], p[(k + 1) % n][1] - p[k][1], p[(k + 1) % n][2] - p[k][2] }, x[3];
        cross3(p[k], d, x);
        ap = fmin(ap, fabs(dot3(x, o->n)) / sqrt(dot3(d, d)));
    }
    o->R = ap - inset;
    double l0 = sqrt(dot3(p[0], p[0]));
    for (int i = 0; i < 3; ++i) o->t[i] = p[0][i] / l0;
    cross3(o->n, o->t, o->b);
    for (int k = 1; k < n; ++k) { double a = atan2(dot3(p[k], o->b), dot3(p[k], o->t)); o->ang[k] = a < 0 ? a + 2 * 3.14159265358979 : a; }
    return o->R > 2 * r;
}

/* the edge k and u along it of angle a */
static void circ_ref_segment(const circ_ref_t* o, double a, int* k, double* u) {
    *k = 0;
    for (int i = 1; i < o->nc; ++i) if (o->ang[i] <= a) *k = i;
    double a0 = o->ang[*k], a1 = *k + 1 < o->nc ? o->ang[*k + 1] : 2 * 3.14159265358979;
    *u = (a - a0) / (a1 - a0);
}

/* within a dash at angle a (by a margin of `ma` radians: > 0 shrinks the dashes, < 0 grows them) */
static bool circ_ref_in_dash(const circ_ref_t* o, double a, double ma) {
    if (!o->dashes) return true;
    for (int k = 0; k < o->nc; ++k) {
        double a0 = o->ang[k], a1 = k + 1 < o->nc ? o->ang[k + 1] : 2 * 3.14159265358979, m = o->dashes;
        for (int j = 0; j < o->dashes; ++j) {
            double s0 = a0 + (a1 - a0) * (j + 0.25) / m, s1 = a0 + (a1 - a0) * (j + 0.75) / m;
            if (a >= s0 + ma && a <= s1 - ma) return true;
        }
    }
    return false;
}

/* the least distance from the ray (from its origin on) to the circle's centre line, over its points in dashes */
static double circ_ref_ray_dist(const circ_ref_t* o, ray_t r, double ma) {
    double best = 1e30, dc;
    sphere_t(r, o->c, 0, &dc);
    if (dc > o->R + 2.0 * o->r + 0.5) return best;          /* nowhere near */
    for (int s = 0; s < 1440; ++s) {
        double a = (s + 0.5) * 2 * 3.14159265358979 / 1440;
        if (!circ_ref_in_dash(o, a, ma)) continue;
        double pt[3];
        for (int i = 0; i < 3; ++i) pt[i] = o->c[i] + o->R * (cos(a) * o->t[i] + sin(a) * o->b[i]);
        double d;
        sphere_t(r, pt, 0, &d);
        double oc[3] = { pt[0] - r.o[0], pt[1] - r.o[1], pt[2] - r.o[2] };
        if (dot3(oc, r.d) < 0) d = sqrt(dot3(oc, oc));
        best = fmin(best, d);
    }
    return best;
}

/* signed distance of x to the circle's tube, and x's angle */
static double circ_ref_sdf(const circ_ref_t* o, const double x[3], double* ang) {
    double d[3] = { x[0] - o->c[0], x[1] - o->c[1], x[2] - o->c[2] };
    double lx = dot3(d, o->t), ly = dot3(d, o->b), lz = dot3(d, o->n);
    double a = atan2(ly, lx);
    *ang = a < 0 ? a + 2 * 3.14159265358979 : a;
    double rho = sqrt(lx * lx + ly * ly) - o->R;
    return sqrt(rho * rho + lz * lz) - o->r;
}

/* POLYGONS with an inset: the circles inside them (md_gl's aromatic rings)
   against tori in double: solid hexagons at several tilts, nearly edge on
   too, and dashed pentagons (the gaps on the corners); every pixel empty or
   covered as the reference says, on the tube (or on a dash's flat end),
   the ids where one is clearly nearest; a solid circle's colour by the
   corners either side (sharpness 1: the nearer one's); a pick */
UTEST(render, polygons_circles_match_their_tori) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 200, H = 150;
    cam_t c = { .tan_half = 0.35f, .aspect = (float)W / H, .n = 0.1f };
    md_render_camera_t cam = make_camera(&c);
    enum { N6 = 4, N5 = 3 };
    double q6[N6][6][3], q5[N5][5][3];
    const double c6[N6][3] = { { -2.6, 1.3, -12 }, { 0.2, 1.4, -12.5 }, { 2.9, 1.2, -11.5 }, { -1.0, -1.2, -12 } };
    const double t6[N6] = { 0, 55, 84, 30 }, s6[N6] = { 0, 20, -10, 100 };
    const double c5[N5][3] = { { 1.8, -1.4, -11.5 }, { -3.2, -1.6, -13 }, { 3.6, -2.2, -13 } };
    const double t5[N5] = { 10, 60, 40 }, s5[N5] = { 0, -40, 70 };
    for (int k = 0; k < N6; ++k) poly_corners(q6[k], 6, c6[k], 1.39, t6[k], s6[k], k == 3, k == 1 ? 0.12 : 0.0);
    for (int k = 0; k < N5; ++k) poly_corners(q5[k], 5, c5[k], 1.18, t5[k], s5[k], false, 0.0);
    float pos6[N6 * 6][3], pos5[N5 * 5][3];
    for (int k = 0; k < N6; ++k) for (int j = 0; j < 6; ++j) for (int i = 0; i < 3; ++i) pos6[6 * k + j][i] = (float)q6[k][j][i];
    for (int k = 0; k < N5; ++k) for (int j = 0; j < 5; ++j) for (int i = 0; i < 3; ++i) pos5[5 * k + j][i] = (float)q5[k][j][i];
    uint32_t col6[N6 * 6];
    const uint32_t pal[6] = { 0xFF0000FFu, 0xFF00FF00u, 0xFFFF0000u, 0xFF00FFFFu, 0xFFFF00FFu, 0xFFFFFF00u };
    for (int i = 0; i < N6 * 6; ++i) col6[i] = pal[i % 6];
    const float R = 0.07f, INS = 0.38f;
    md_render_batch_t b[2];
    memset(b, 0, sizeof(b));
    b[0].key = 61; b[0].kind = MD_RENDER_POLYGONS; b[0].count = N6; b[0].points.count = N6 * 6;
    b[0].points.position = f_buf(&fx, pos6, sizeof(pos6));
    b[0].points.color    = f_buf(&fx, col6, sizeof(col6));
    b[0].color.mode = MD_RENDER_COLOR_POINTS; b[0].color.sharpness = 1.0f;
    b[0].style.polygons.corners = 6; b[0].style.polygons.radius = R; b[0].style.polygons.inset = INS;
    b[1].key = 62; b[1].kind = MD_RENDER_POLYGONS; b[1].count = N5; b[1].points.count = N5 * 5;
    b[1].points.position = f_buf(&fx, pos5, sizeof(pos5));
    b[1].style.polygons.corners = 5; b[1].style.polygons.radius = R; b[1].style.polygons.inset = INS; b[1].style.polygons.dashes = 2;
    md_render_frame_t f = frame_of(W, H, cam, b, 2);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));

    circ_ref_t ref[N6 + N5];
    int nref = 0;
    for (int k = 0; k < N6; ++k) { ASSERT_TRUE(circ_ref_build(&ref[nref], (const double (*)[3])q6[k], 6, INS, R, 0)); ref[nref].key = 61; ref[nref].linear = (uint32_t)k; nref++; }
    for (int k = 0; k < N5; ++k) { ASSERT_TRUE(circ_ref_build(&ref[nref], (const double (*)[3])q5[k], 5, INS, R, 2)); ref[nref].key = 62; ref[nref].linear = (uint32_t)k; nref++; }

    md_render_targets_t tg = md_render_view_targets(fx.view);
    uint32_t* id  = malloc(W * H * 8);
    float*    dp  = malloc(W * H * 4);
    uint32_t* col = malloc(W * H * 4);
    ASSERT_TRUE(f_read(&fx, tg.id, id, W * H * 8) && f_read(&fx, tg.depth, dp, W * H * 4) && f_read(&fx, tg.color, col, W * H * 4));
    int bad = 0, covered = 0, judged = 0, ends = 0, coloured = 0, gaps_seen = 0;
    double worst = 0;
    for (uint32_t y = 0; y < H; ++y) {
        for (uint32_t x = 0; x < W; ++x) {
            const uint32_t i = y * W + x;
            ray_t r = pixel_ray(&c, W, H, x + 0.5, y + 0.5);
            bool hit_in = false, hit_out = false;
            for (int k = 0; k < nref; ++k) {
                const double px = 12.0 * 2.0 * c.tan_half / H;
                const double grow = (ref[k].r + px) / (ref[k].R - ref[k].r);   /* past a dash's end: its tube, a pixel */
                hit_in  |= circ_ref_ray_dist(&ref[k], r, ref[k].dashes ? 1.5 * grow : 0.0) < ref[k].r - 0.6 * px;
                hit_out |= circ_ref_ray_dist(&ref[k], r, ref[k].dashes ? -grow : 0.0) < ref[k].r + 0.6 * px;
            }
            if (id[2 * i + 1] == 0) {
                if (hit_in) { bad++; if (bad < 8) printf("    pixel (%u, %u): empty, a circle there\n", x, y); }
                continue;
            }
            covered++;
            if (!hit_out) { bad++; if (bad < 8) printf("    pixel (%u, %u): covered, no circle near\n", x, y); continue; }
            const double zv = -c.n / dp[i], tp = zv / r.d[2];
            const double P[3] = { r.o[0] + r.d[0] * tp, r.o[1] + r.d[1] * tp, r.o[2] + r.d[2] * tp };
            const double pxw = -P[2] * 2.0 * c.tan_half / H;
            int who = -1;
            double best = 1e30, second = 1e30, ang = 0;
            for (int k = 0; k < nref; ++k) {
                double a, sd = circ_ref_sdf(&ref[k], P, &a);
                if (ref[k].dashes && !circ_ref_in_dash(&ref[k], a, -0.02)) sd = fabs(sd) + 1.0;     /* in a gap: not this one */
                if (fabs(sd) < fabs(best)) { second = best; best = sd; who = k; ang = a; } else if (fabs(sd) < fabs(second)) second = sd;
            }
            bool on_end = false;
            if (who >= 0 && ref[who].dashes && best < -0.05 * pxw) {
                /* inside the tube: on a dash's flat end, at the end's angle */
                on_end = !circ_ref_in_dash(&ref[who], ang, 0.02) && circ_ref_in_dash(&ref[who], ang, -0.02);
                ends += on_end;
            }
            worst = fmax(worst, on_end ? 0.0 : fabs(best) / pxw);
            if (!on_end && fabs(best) > 0.05 * pxw) { bad++; if (bad < 8) printf("    pixel (%u, %u): %.3f px off the circles\n", x, y, best / pxw); continue; }
            if (fabs(second) - fabs(best) > pxw) {
                judged++;
                const uint32_t slot = (id[2 * i + 1] >> 20) - 1u;
                if (id[2 * i] != ref[who].linear || slot != ref[who].key - 61u) {
                    bad++;
                    if (bad < 8) printf("    pixel (%u, %u): batch %u element %u, expected %u %u\n", x, y, slot, id[2 * i], ref[who].key - 61u, ref[who].linear);
                    continue;
                }
                if (!ref[who].dashes) {
                    int k; double u;
                    circ_ref_segment(&ref[who], ang, &k, &u);
                    if (fabs(u - 0.5) > 0.06) {
                        const int kk = u < 0.5 ? k : (k + 1) % ref[who].nc;
                        coloured++;
                        if (!near_rgba(col[i] | 0xFF000000u, pal[kk])) { bad++; if (bad < 8) printf("    pixel (%u, %u): colour %08x, expected %08x\n", x, y, col[i], pal[kk]); }
                    }
                }
            }
        }
    }
    /* the gaps: along a dashed circle's centre line, a gap's middle is empty */
    {
        const circ_ref_t* d = &ref[N6];
        for (int k = 0; k < d->nc; ++k) {
            const double a = d->ang[k] + (k ? 0.0 : 0.0);        /* a corner: a gap's middle */
            double pt[3];
            for (int i = 0; i < 3; ++i) pt[i] = d->c[i] + d->R * (cos(a) * d->t[i] + sin(a) * d->b[i]);
            const uint32_t qx = (uint32_t)((pt[0] / -pt[2] / (c.tan_half * c.aspect) * 0.5 + 0.5) * W);
            const uint32_t qy = (uint32_t)((0.5 - pt[1] / -pt[2] / c.tan_half * 0.5) * H);
            gaps_seen += id[2 * (qy * W + qx) + 1] == 0;
        }
    }
    printf("    %d pixels covered, %d of one circle clearly (%d on dashes' ends, %d coloured); the farthest %.4f px off; %d of 5 gaps empty\n",
           covered, judged, ends, coloured, worst, gaps_seen);
    EXPECT_EQ(0, bad);
    EXPECT_GT(covered, 800);
    EXPECT_GT(judged, covered / 2);
    EXPECT_GT(coloured, 300);
    if (getenv("MD_RENDER_DUMP")) {
        FILE* fp = fopen(getenv("MD_RENDER_DUMP"), "wb");
        if (fp) {
            fprintf(fp, "P6\n%u %u\n255\n", W, H);
            for (uint32_t i = 0; i < W * H; ++i) { uint8_t px3[3] = { (uint8_t)col[i], (uint8_t)(col[i] >> 8), (uint8_t)(col[i] >> 16) }; fwrite(px3, 1, 3, fp); }
            fclose(fp);
        }
    }
    EXPECT_EQ(5, gaps_seen);
    /* a pick on the first circle, at its first corner's angle, past the gap: that circle */
    {
        const circ_ref_t* d = &ref[0];
        const double a = 0.5 * d->ang[1];
        double pt[3];
        for (int i = 0; i < 3; ++i) pt[i] = d->c[i] + d->R * (cos(a) * d->t[i] + sin(a) * d->b[i]) + d->n[i] * 0.5 * d->r;
        const uint32_t qx = (uint32_t)((pt[0] / -pt[2] / (c.tan_half * c.aspect) * 0.5 + 0.5) * W);
        const uint32_t qy = (uint32_t)((0.5 - pt[1] / -pt[2] / c.tan_half * 0.5) * H);
        md_render_pick_t pk;
        ASSERT_TRUE(f_pick(&fx, qx, qy, &pk));
        EXPECT_EQ(61u, (unsigned)pk.key);
        EXPECT_EQ(0u, pk.element);
    }
    free(id); free(dp); free(col);
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
    b.points.marks = f_marks(&fx, &sel, NULL, 32);
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

/* The same tiers for a batch drawn as instances (periodic images,
   assemblies): an unchanged frame draws nothing, a change of state or
   appearance no raster; their transforms and their count are geometry. (Every
   draw of an instanced batch was a raster: the view keeps a batch's
   transforms apart from it, and the comparison read their count off the
   pointer it no longer had.) */
UTEST(render, lazy_tiers_with_instances) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 32, H = 32;
    cam_t c = { .tan_half = 0.5f, .aspect = 1.0f, .n = 0.1f };
    float pos[3] = { 0, 0, -5 };
    uint32_t sel = 1;
    md_render_batch_t b = {0};
    b.key = 1; b.kind = MD_RENDER_SPHERES; b.count = 1; b.points.count = 1;
    b.points.position = f_buf(&fx, pos, sizeof(pos));
    b.points.marks = f_marks(&fx, &sel, NULL, 32);
    b.style.spheres.radius_add = 0.5f;
    b.color.mode = MD_RENDER_COLOR_UNIFORM;
    b.color.uniform = 0xFFFFFFFFu;
    float xf[32] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, -1,0,0,1,   1,0,0,0, 0,1,0,0, 0,0,1,0, 1,0,0,1 };
    b.transform = xf;
    b.transform_count = 2;
    md_render_frame_t f = frame_of(W, H, make_camera(&c), &b, 1);

    md_render_status_t st;
    const md_render_work_t ALL = MD_RENDER_WORK_RASTER | MD_RENDER_WORK_RESOLVE | MD_RENDER_WORK_MASK;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    EXPECT_EQ(ALL, st.work);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    EXPECT_EQ(0u, st.work);                                          /* unchanged */
    f.content[MD_RENDER_CONTENT_STATE]++;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    EXPECT_EQ(MD_RENDER_WORK_MASK, st.work);
    f.content[MD_RENDER_CONTENT_APPEARANCE]++;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    EXPECT_EQ(MD_RENDER_WORK_RESOLVE | MD_RENDER_WORK_MASK, st.work);
    b.color.uniform = 0xFF0000FFu;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    EXPECT_EQ(MD_RENDER_WORK_RESOLVE | MD_RENDER_WORK_MASK, st.work);

    /* A transform moved: a raster; then nothing again */
    xf[28] = 1.5f;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    EXPECT_EQ(ALL, st.work);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    EXPECT_EQ(0u, st.work);
    /* One instance fewer, and one more */
    b.transform_count = 1;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    EXPECT_EQ(ALL, st.work);
    b.transform_count = 2;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    EXPECT_EQ(ALL, st.work);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    EXPECT_EQ(0u, st.work);
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

/* The motion of a pixel is that of the surface that makes the most of it: a
   moving transparent sphere in front of a still opaque one has its own motion
   where it is nearly opaque (0.8 of the pixel against 0.2 behind it), and
   leaves the opaque one's (0) where it is faint (0.2 against 0.8) */
static float half_to_float(uint32_t h);

UTEST(render, motion_of_the_dominant_surface) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    md_gpu_device_info_t info;
    md_gpu_device_info(fx.dev, &info);
    if (!info.fragment_stores) UTEST_SKIP("no stores in fragment stages (the K-buffer's)");
    const uint32_t W = 64, H = 64;
    cam_t c = { .tan_half = 0.5f, .aspect = 1.0f, .n = 0.1f };
    md_render_camera_t cam = make_camera(&c);
    const float vel[3] = { 2.0f, -1.0f, 0.0f };
    const double dt = 0.25;
    uint32_t pv[2];
    md_render_pack_velocity(vel, pv);
    const float opq[3] = { 0.0f, 0.0f, -20.0f };
    const float tp0[3] = { 0.0f, 0.0f, -10.0f };
    const float tp1[3] = { (float)(tp0[0] + vel[0] * dt), (float)(tp0[1] + vel[1] * dt), tp0[2] };
    /* the moved sphere's front point (depth 8.5) under its centre's pixel, and its motion there */
    const uint32_t px = (uint32_t)(tp1[0] / 10.0 / c.tan_half * 0.5 * W + W / 2);
    const uint32_t py = (uint32_t)(H / 2 - tp1[1] / 10.0 / c.tan_half * 0.5 * H);
    const double ex = vel[0] * dt / 8.5 / c.tan_half * 0.5 * W;
    const double ey = -vel[1] * dt / 8.5 / c.tan_half * 0.5 * H;

    md_gpu_addr_t O = f_buf(&fx, opq, sizeof(opq));
    md_gpu_addr_t T = f_buf(&fx, tp0, sizeof(tp0));
    md_render_batch_t b[2];
    memset(b, 0, sizeof(b));
    b[0].key = 1; b[0].kind = MD_RENDER_SPHERES; b[0].count = 1; b[0].points.count = 1; b[0].points.position = O;
    b[0].style.spheres.radius_add = 6.0f;
    b[1].key = 2; b[1].kind = MD_RENDER_SPHERES; b[1].count = 1; b[1].points.count = 1; b[1].points.position = T;
    b[1].points.velocity = f_buf(&fx, pv, sizeof(pv));
    b[1].style.spheres.radius_add = 1.5f;
    uint32_t* mo = malloc(W * H * 4);
    static const float transparency[2] = { 0.2f, 0.8f };
    for (int i = 0; i < 2; ++i) {
        b[1].color.transparency = transparency[i];
        md_gpu_upload(fx.gfx, T, tp0, sizeof(tp0));
        md_render_frame_t f = frame_of(W, H, cam, b, 2);
        f.time = 0.0;
        ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
        md_gpu_upload(fx.gfx, T, tp1, sizeof(tp1));
        f.content[MD_RENDER_CONTENT_GEOMETRY]++;
        f.time = dt;
        ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
        ASSERT_TRUE(f_read(&fx, md_render_view_targets(fx.view).motion, mo, W * H * 4));
        const uint32_t w = mo[py * W + px];
        const double mx = half_to_float(w & 0xFFFFu), my = half_to_float(w >> 16);
        if (i == 0) {
            EXPECT_NEAR(ex, mx, 0.15 * fabs(ex));
            EXPECT_NEAR(ey, my, 0.15 * fabs(ey));
        } else {
            EXPECT_EQ(0.0, mx);
            EXPECT_EQ(0.0, my);
        }
        /* outside the transparent sphere: the opaque one's */
        EXPECT_EQ(0u, mo[(H / 2 + 12) * W + W / 2 - 12]);
    }
    free(mo);
    f_close(&fx);
}

/* The marks' tags: a set mask channel holds (8 + tag) / 15, the visible
   surface's tag in the visible channels and the largest in the anywhere
   ones; a capsule has a state when both its points have it, with the larger
   tag. The export carries the tag of state 1 in the colour's alpha. And the
   CPU helpers that write marks from sets. */
UTEST(render, masks_carry_tags) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 64, H = 32;
    cam_t c = { .ortho = true, .half_h = 2.0f, .aspect = 2.0f, .n = 0.0f, .far_ = 100.0f };
    /* 0 highlighted (tag 2) in front of 1 highlighted (tag 6); 2 selected (tag 5);
       capsules 3-4 highlighted (tags 1, 4) and 5-6 (only 5 highlighted) */
    float pos[21] = { -1.0f, 0, -10,   0.0f, 0, -20,   2.5f, 0, -10,
                      -3.0f, 1.5f, -10,  -1.0f, 1.5f, -10,   1.0f, 1.5f, -10,   3.0f, 1.5f, -10 };
    uint8_t marks[7] = {
        MD_RENDER_MARK_STATE1 | MD_RENDER_MARK_TAG(2), MD_RENDER_MARK_STATE1 | MD_RENDER_MARK_TAG(6),
        MD_RENDER_MARK_STATE0 | MD_RENDER_MARK_TAG(5),
        MD_RENDER_MARK_STATE1 | MD_RENDER_MARK_TAG(1), MD_RENDER_MARK_STATE1 | MD_RENDER_MARK_TAG(4),
        MD_RENDER_MARK_STATE1 | MD_RENDER_MARK_TAG(3), MD_RENDER_MARK_TAG(7) };
    uint32_t index[4] = { 3, 4, 5, 6 };
    md_gpu_addr_t P = f_buf(&fx, pos, sizeof(pos)), M = f_buf(&fx, marks, sizeof(marks));
    md_render_batch_t b[2];
    memset(b, 0, sizeof(b));
    b[0].key = 1; b[0].kind = MD_RENDER_SPHERES; b[0].count = 3;
    b[0].style.spheres.radius_add = 0.9f;
    b[1].key = 2; b[1].kind = MD_RENDER_CAPSULES; b[1].count = 2;
    b[1].index = f_buf(&fx, index, sizeof(index));
    b[1].style.capsules.radius = 0.3f;
    for (int k = 0; k < 2; ++k) { b[k].points.position = P; b[k].points.count = 7; b[k].points.marks = M; }
    md_render_frame_t f = frame_of(W, H, make_camera(&c), b, 2);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    uint32_t* m = malloc(W * H * 4);
    ASSERT_TRUE(f_read(&fx, md_render_view_targets(fx.view).mask, m, W * H * 4));
    #define AT(wx, wy) m[(uint32_t)((0.5 - (wy) / 4.0) * H) * W + (uint32_t)(((wx) / 4.0 * 0.5 + 0.5) * W)]
    #define CH(tag) ((8u + (tag)) * 17u)
    EXPECT_EQ((CH(6) << 24) | (CH(2) << 16), AT(-0.5, 0.0));   /* 0 in view (tag 2), 1 behind it: the larger, 6 */
    EXPECT_EQ((CH(6) << 24) | (CH(6) << 16), AT(0.5, 0.0));    /* 1 in view */
    EXPECT_EQ((CH(5) << 8) | CH(5), AT(2.5, 0.0));             /* the selection's tag */
    EXPECT_EQ((CH(4) << 24) | (CH(4) << 16), AT(-2.0, 1.5));   /* both points highlighted: the larger tag */
    EXPECT_EQ(0u, AT(2.0, 1.5));                               /* one of them: not */

    /* The export: the state 1 tag of the visible surface, else of the largest */
    md_gpu_texture_desc_t td = {0};
    td.type = MD_GPU_TEX_2D; td.usage = MD_GPU_TEX_STORAGE; td.width = W; td.height = H;
    td.format = MD_GPU_FORMAT_RGBA8_UNORM; td.label = "export color";
    md_gpu_texture_t xc = md_gpu_texture_create(fx.gfx, &td);
    md_render_export_t ex = { 0, xc, 0 };
    ASSERT_TRUE(md_render_view_export(fx.gfx, fx.view, &ex));
    ASSERT_TRUE(f_read(&fx, xc, m, W * H * 4));
    EXPECT_EQ(12u | (2u << 4) | 128u, AT(-0.5, 0.0) >> 24);        /* 128: of a batch with marks */
    EXPECT_EQ(12u | (6u << 4) | 128u, AT(0.5, 0.0) >> 24);
    EXPECT_EQ(3u | 128u, AT(2.5, 0.0) >> 24);
    #undef CH
    #undef AT
    md_gpu_texture_destroy(xc);
    free(m);

    /* The helpers: states from sets, then tags (the later wins) */
    uint64_t w0 = 0x5ull, w1 = 0x6ull, wt = 0xCull, wu = 0x8ull;
    md_bitset_t s0 = { .words = &w0, .beg = 0, .end = 64 }, s1 = { .words = &w1, .beg = 0, .end = 64 };
    md_bitset_t st = { .words = &wt, .beg = 0, .end = 64 }, su = { .words = &wu, .beg = 0, .end = 64 }, none = {0};
    uint8_t out[5];
    memset(out, 0xFF, sizeof(out));
    md_render_marks_write(out, 4, s0, s1);
    EXPECT_EQ(MD_RENDER_MARK_STATE0, out[0]);
    EXPECT_EQ(MD_RENDER_MARK_STATE1, out[1]);
    EXPECT_EQ(MD_RENDER_MARK_STATE0 | MD_RENDER_MARK_STATE1, out[2]);
    EXPECT_EQ(0, out[3]);
    EXPECT_EQ(0xFF, out[4]);                                   /* beyond the count: untouched */
    md_render_marks_tag(out, 4, st, 3);
    md_render_marks_tag(out, 4, su, 6);
    EXPECT_EQ(MD_RENDER_MARK_STATE0 | MD_RENDER_MARK_STATE1 | MD_RENDER_MARK_TAG(3), out[2]);
    EXPECT_EQ(MD_RENDER_MARK_TAG(6), out[3]);
    md_render_marks_write(out, 4, none, none);
    EXPECT_EQ(0, out[0] | out[1] | out[2] | out[3]);
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
    b.points.marks = f_marks(&fx, &sel, &hl, 32);
    b.style.spheres.radius_add = 0.9f;
    md_render_frame_t f = frame_of(W, H, make_camera(&c), &b, 1);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    uint32_t* m = malloc(W * H * 4);
    ASSERT_TRUE(f_read(&fx, md_render_view_targets(fx.view).mask, m, W * H * 4));
    #define AT(wx) m[(H / 2) * W + (uint32_t)(((wx) / 4.0 * 0.5 + 0.5) * W)]
    EXPECT_EQ(0x00008800u | 0x00000088u, AT(-1.5));   /* selection, visible and there   */
    EXPECT_EQ(0x00008800u, AT(-0.3));                 /* selection hidden behind point 1 */
    EXPECT_EQ(0x88000000u | 0x00880000u, AT(2.5));    /* highlight, visible             */
    EXPECT_EQ(0u, AT(0.6));                           /* point 1: nothing               */
    #undef AT

    /* Selection moves to point 1: STATE only */
    sel = 1u << 1;
    f_marks_upload(&fx, b.points.marks, &sel, &hl, 32);
    f.content[MD_RENDER_CONTENT_STATE]++;
    md_render_status_t st;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    EXPECT_EQ(MD_RENDER_WORK_MASK, st.work);
    ASSERT_TRUE(f_read(&fx, md_render_view_targets(fx.view).mask, m, W * H * 4));
    #define AT(wx) m[(H / 2) * W + (uint32_t)(((wx) / 4.0 * 0.5 + 0.5) * W)]
    EXPECT_EQ(0u, AT(-1.5));
    EXPECT_EQ(0x00008800u | 0x00000088u, AT(-0.3));
    #undef AT

    /* The camera moves a little: a raster that keeps the state lists, and the
       masks still show the selection */
    f.camera.view[12] = 0.05f;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
    EXPECT_EQ(MD_RENDER_WORK_RASTER | MD_RENDER_WORK_RESOLVE | MD_RENDER_WORK_MASK, st.work);
    ASSERT_TRUE(f_read(&fx, md_render_view_targets(fx.view).mask, m, W * H * 4));
    #define AT(wx) m[(H / 2) * W + (uint32_t)(((wx) / 4.0 * 0.5 + 0.5) * W)]
    EXPECT_EQ(0u, AT(-1.5));
    EXPECT_EQ(0x00008800u | 0x00000088u, AT(-0.3));
    EXPECT_EQ(0x88000000u | 0x00880000u, AT(2.5));    /* the highlight too */
    #undef AT
    free(m);
    f_close(&fx);
}

/* The export copies the targets into the caller's textures, bit for bit:
   depth, the colour with the mask's channels as bits in its alpha, the normal
   and the motion together */
UTEST(render, export_copies_the_targets) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 64, H = 32;
    cam_t c = { .tan_half = 0.4f, .aspect = 2.0f, .n = 0.1f };
    /* Point 0 selected, behind point 1 which overlaps it; point 2 highlighted and moving */
    float pos[9] = { -1.0f, 0, -20,   0.0f, 0, -12,   2.5f, 0, -10 };
    float vel[9] = { 0 }; vel[6] = 3.0f; vel[7] = 1.0f;
    uint32_t pv[6];
    for (int i = 0; i < 3; ++i) md_render_pack_velocity(vel + 3 * i, pv + 2 * i);
    uint32_t sel = 1u << 0, hl = 1u << 2;
    uint32_t col[3] = { 0xFF2040E0u, 0xFF10C010u, 0xFFE0E0E0u };
    md_render_batch_t b = {0};
    b.key = 1; b.kind = MD_RENDER_SPHERES; b.count = 3; b.points.count = 3;
    b.points.position = f_buf(&fx, pos, sizeof(pos));
    b.points.velocity = f_buf(&fx, pv, sizeof(pv));
    b.points.color    = f_buf(&fx, col, sizeof(col));
    b.points.marks = f_marks(&fx, &sel, &hl, 32);
    b.style.spheres.radius_add = 1.2f;
    /* and a sphere of a batch without marks: none of the states, nor the bit
       of a batch with marks */
    float pos2[3] = { -4.0f, 0, -12 };
    md_render_batch_t bb[2] = { b, b };
    bb[1].key = 2; bb[1].count = 1; bb[1].points.count = 1;
    bb[1].points.position = f_buf(&fx, pos2, sizeof(pos2));
    bb[1].points.velocity = 0;
    bb[1].points.marks = 0;
    md_render_frame_t f = frame_of(W, H, make_camera(&c), bb, 2);
    f.time = 0.0;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    f.time = 0.1;
    f.camera.view[12] = 0.1f;                  /* and the camera moves: motion everywhere */
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));

    md_gpu_texture_desc_t td = {0};
    td.type = MD_GPU_TEX_2D; td.usage = MD_GPU_TEX_STORAGE; td.width = W; td.height = H;
    td.format = MD_GPU_FORMAT_R32_FLOAT;    td.label = "export depth";   md_gpu_texture_t xd = md_gpu_texture_create(fx.gfx, &td);
    td.format = MD_GPU_FORMAT_RGBA8_UNORM;  td.label = "export color";   md_gpu_texture_t xc = md_gpu_texture_create(fx.gfx, &td);
    td.format = MD_GPU_FORMAT_RGBA16_FLOAT; td.label = "export surface"; md_gpu_texture_t xs = md_gpu_texture_create(fx.gfx, &td);
    ASSERT_TRUE(xd && xc && xs);
    md_render_export_t ex = { xd, xc, xs };
    ASSERT_TRUE(md_render_view_export(fx.gfx, fx.view, &ex));

    /* A texture of another size or format is refused */
    md_render_export_t wrong = { xc, 0, 0 };
    EXPECT_FALSE(md_render_view_export(fx.gfx, fx.view, &wrong));

    md_render_targets_t t = md_render_view_targets(fx.view);
    const size_t n = (size_t)W * H;
    uint32_t* d0 = malloc(n * 4); uint32_t* d1 = malloc(n * 4);
    uint32_t* c0 = malloc(n * 4); uint32_t* c1 = malloc(n * 4); uint32_t* m0 = malloc(n * 4);
    uint32_t* n0 = malloc(n * 4); uint32_t* v0 = malloc(n * 4); uint32_t* s1 = malloc(n * 8);
    uint32_t* id = malloc(n * 8);
    ASSERT_TRUE(f_read(&fx, t.id, id, n * 8));
    ASSERT_TRUE(f_read(&fx, t.depth, d0, n * 4) && f_read(&fx, xd, d1, n * 4));
    ASSERT_TRUE(f_read(&fx, t.color, c0, n * 4) && f_read(&fx, xc, c1, n * 4) && f_read(&fx, t.mask, m0, n * 4));
    ASSERT_TRUE(f_read(&fx, t.normal, n0, n * 4) && f_read(&fx, t.motion, v0, n * 4) && f_read(&fx, xs, s1, n * 8));
    int bad_depth = 0, bad_color = 0, bad_mask = 0, bad_surface = 0, covered = 0, masked = 0, moving = 0, unmarked = 0;
    for (size_t i = 0; i < n; ++i) {
        bad_depth += d0[i] != d1[i];
        bad_color += (c0[i] & 0xFFFFFFu) != (c1[i] & 0xFFFFFFu);
        uint32_t bits = 0;
        for (int k = 0; k < 4; ++k) bits |= ((m0[i] >> (8 * k)) & 0xFFu) >= 128u ? 1u << k : 0u;
        const uint32_t slot = id[2 * i + 1] >> 20;     /* the batch + 1 */
        if (slot == 1) bits |= 128u;
        unmarked += slot == 2;
        bad_mask += (c1[i] >> 24) != bits;
        bad_surface += s1[2 * i] != n0[i] || s1[2 * i + 1] != v0[i];
        covered += d0[i] != 0;
        masked  += bits != 0;
        moving  += v0[i] != 0;
    }
    EXPECT_EQ(0, bad_depth);
    EXPECT_EQ(0, bad_color);
    EXPECT_EQ(0, bad_mask);
    EXPECT_EQ(0, bad_surface);
    EXPECT_GT(covered, 100);
    EXPECT_GT(masked, 20);
    EXPECT_GT(moving, 20);
    EXPECT_GT(unmarked, 20);
    free(d0); free(d1); free(c0); free(c1); free(m0); free(n0); free(v0); free(s1); free(id);
    md_gpu_texture_destroy(xd); md_gpu_texture_destroy(xc); md_gpu_texture_destroy(xs);
    f_close(&fx);
}

/* A jitter the caller gives: the raster is that of the camera moved by it (as
   md_render moves it for its own), the motion is without it, and the view
   rasterises again when, and only when, it changes */
UTEST(render, jitter_given) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 48, H = 40;
    cam_t c = { .tan_half = 0.4f, .aspect = (float)W / H, .n = 0.1f };
    md_render_camera_t cam = make_camera(&c);
    const double ctr[2][3] = { { 0.3, 0.1, -10 }, { -1.2, -0.4, -7 } };
    const double rad[2]    = { 1.6, 0.9 };
    float pos[6], radius[2];
    for (int i = 0; i < 2; ++i) { for (int k = 0; k < 3; ++k) pos[3 * i + k] = (float)ctr[i][k]; radius[i] = (float)rad[i]; }
    md_render_batch_t b = {0};
    b.key = 3; b.kind = MD_RENDER_SPHERES; b.count = 2; b.points.count = 2;
    b.points.position = f_buf(&fx, pos, sizeof(pos));
    b.points.radius = f_buf(&fx, radius, sizeof(radius));
    b.style.spheres.radius_scale = 1.0f;

    md_render_view_t ref = md_render_view_create(fx.r, "jitter reference");
    ASSERT_TRUE(ref != NULL);
    const float off[2][2] = { { 0.31f, -0.27f }, { -0.42f, 0.18f } };
    for (int j = 0; j < 2; ++j) {
        md_render_frame_t f = frame_of(W, H, cam, &b, 1);
        f.jitter = MD_RENDER_JITTER_GIVEN;
        f.jitter_offset[0] = off[j][0];
        f.jitter_offset[1] = off[j][1];
        md_render_status_t st;
        ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
        EXPECT_TRUE((st.work & MD_RENDER_WORK_RASTER) != 0);
        EXPECT_EQ(off[j][0], st.jitter[0]);
        EXPECT_EQ(off[j][1], st.jitter[1]);
        /* The same offset again: nothing to do */
        ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, &st));
        EXPECT_EQ(0u, st.work);

        /* The reference: the camera moved as md_render moves it, unjittered */
        md_render_frame_t g = frame_of(W, H, cam, &b, 1);
        const float ndx = 2.0f * off[j][0] / (float)W, ndy = -2.0f * off[j][1] / (float)H;
        for (int col = 0; col < 4; ++col) {
            g.camera.proj[col * 4 + 0] += ndx * g.camera.proj[col * 4 + 3];
            g.camera.proj[col * 4 + 1] += ndy * g.camera.proj[col * 4 + 3];
        }
        ASSERT_TRUE(md_render_view_draw(fx.gfx, ref, &g, NULL));
        const size_t n = (size_t)W * H;
        uint32_t* i0 = malloc(n * 8); uint32_t* i1 = malloc(n * 8);
        float* p0 = malloc(n * 4); float* p1 = malloc(n * 4);
        ASSERT_TRUE(f_read(&fx, md_render_view_targets(fx.view).id, i0, n * 8) && f_read(&fx, md_render_view_targets(ref).id, i1, n * 8));
        ASSERT_TRUE(f_read(&fx, md_render_view_targets(fx.view).depth, p0, n * 4) && f_read(&fx, md_render_view_targets(ref).depth, p1, n * 4));
        int bad = 0, covered = 0;
        for (size_t i = 0; i < n; ++i) {
            bad += i0[2 * i] != i1[2 * i] || (i0[2 * i + 1] >> 20 != 0) != (i1[2 * i + 1] >> 20 != 0) || p0[i] != p1[i];
            covered += (i0[2 * i + 1] >> 20) != 0;
        }
        EXPECT_EQ(0, bad);
        EXPECT_GT(covered, (int)(n / 10));
        free(i0); free(i1); free(p0); free(p1);
    }
    /* Still geometry and camera: the motion of the second raster is without
       the jitter's change, nothing (half-float 0, or -0 from rounding) */
    uint32_t* mo = malloc((size_t)W * H * 4);
    ASSERT_TRUE(f_read(&fx, md_render_view_targets(fx.view).motion, mo, (size_t)W * H * 4));
    float worst = 0.0f;
    for (size_t i = 0; i < (size_t)W * H; ++i) {
        for (int k = 0; k < 2; ++k) {
            uint32_t h = (mo[i] >> (16 * k)) & 0xFFFFu, sgn = h >> 15, e = (h >> 10) & 31u, m = h & 1023u, bits;
            bits = e ? ((sgn << 31) | ((e + 112u) << 23) | (m << 13)) : (sgn << 31);
            float v; memcpy(&v, &bits, 4);
            worst = fmaxf(worst, fabsf(v));
        }
    }
    EXPECT_LT(worst, 1e-3f);
    free(mo);
    md_render_view_destroy(fx.gfx, ref);
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
    md_render_encode_id(1, 3 * 1000 + 17, id);
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

/* md_render's occlusion settings, for a test of the path its defaults skip,
   and its own choice whether to cull by occlusion at all (in md_render.c, not
   in its header) */

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
    md_render_debug_occlusion_auto(false);      /* the culling itself, every raster: not the view's choice of it */
    const uint32_t W = 160, H = 120;
    cam_t c = { .tan_half = 0.6f, .aspect = (float)W / H, .n = 0.1f };
    dense_t d;
    dense_scene(&fx, &d);

    md_render_view_t ref = md_render_view_create(fx.r, "reference");
    md_render_frame_t fr = frame_of(W, H, make_camera(&c), d.b, 2);
    md_render_debug_view(ref, MD_RENDER_DEBUG_NO_OCCLUSION | MD_RENDER_DEBUG_NO_COMPUTE_RASTER, 0.0f);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, ref, &fr, NULL));

    md_render_frame_t f = frame_of(W, H, make_camera(&c), d.b, 2);
    md_render_debug_view(fx.view, 0, 4.0f);
    md_render_stats_t st;
    /* First raster: no pyramid yet, so nothing is occluded in phase 1 */
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    md_gpu_stream_sync(fx.gfx);
    ASSERT_TRUE(md_render_view_stats(fx.view, &st));
    EXPECT_EQ((uint64_t)(d.n_sph + d.n_cap), st.elements);
    EXPECT_EQ(0u, st.occlusion_retested);
    if (info.atomic_int64_min_max) EXPECT_GT(st.compute[0], 0u);
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
        md_render_debug_view(vo, 0, 4.0f);
        ASSERT_TRUE(md_render_view_draw(fx.gfx, vo, &fo, NULL));
        fo.content[MD_RENDER_CONTENT_GEOMETRY]++;
        ASSERT_TRUE(md_render_view_draw(fx.gfx, vo, &fo, NULL));
        md_gpu_stream_sync(fx.gfx);
        md_render_stats_t so;
        ASSERT_TRUE(md_render_view_stats(vo, &so));
        EXPECT_GT(so.occluded, (d.n_sph + d.n_cap) / 4);
        md_render_view_t ro = md_render_view_create(fx.r, "its reference");
        md_render_frame_t fro = fo;
        md_render_debug_view(ro, MD_RENDER_DEBUG_NO_OCCLUSION | MD_RENDER_DEBUG_NO_COMPUTE_RASTER, 0.0f);
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
        md_render_debug_view(vj, 0, 4.0f);
        md_render_frame_t fjr = fj;
        md_render_debug_view(rj, MD_RENDER_DEBUG_NO_OCCLUSION | MD_RENDER_DEBUG_NO_COMPUTE_RASTER, 0.0f);
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
    md_gpu_addr_t S = f_buf(&fx, NULL, words * 32);
    md_gpu_memset(fx.gfx, S, 0, words * 32);
    d.b[0].points.marks = S;
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));   /* the address is geometry: a raster */
    f_marks_upload(&fx, S, sel, NULL, words * 32);
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
    md_render_debug_occlusion_auto(true);
    f_close(&fx);
}

/* A view culls by occlusion only while it pays (occlusion_wanted): a raster
   that tested it and culled little turns it off for the rasters after,
   until the elements drawn more than double. Whether a raster culled is
   whether it built a pyramid ("hiz1"). Nothing of the image changes. */
static void count_hiz(void* user, const char* pass) { if (!strcmp(pass, "hiz1")) ++*(int*)user; }

static int raster_with_hiz(fix_t* fx, md_render_view_t v, md_render_frame_t* f, int* hiz) {
    const int before = *hiz;
    f->content[MD_RENDER_CONTENT_GEOMETRY]++;
    if (!md_render_view_draw(fx->gfx, v, f, NULL)) return -1;
    md_gpu_stream_sync(fx->gfx);        /* its counts complete: the next raster reads them */
    return *hiz - before;
}

UTEST(render, occlusion_off_where_it_culls_little) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 160, H = 120;
    cam_t c = { .tan_half = 0.6f, .aspect = (float)W / H, .n = 0.1f };
    /* sparse: 10 x 8 spheres apart in one plane, nothing behind anything;
       dense: a block of 32 x 24 x 24 overlapping ones, the front layer
       hiding the rest. All in the hardware raster (no compute raster), so
       every culled element weighs fully. */
    enum { SX = 10, SY = 8, DX = 32, DY = 24, DZ = 24, NS = SX * SY, ND = DX * DY * DZ };
    float* sp = malloc(NS * 12);
    float* dp = malloc((size_t)ND * 12);
    for (int y = 0; y < SY; ++y) for (int x = 0; x < SX; ++x) {
        float* q = sp + 3 * (y * SX + x);
        q[0] = (x - SX / 2) * 3.0f; q[1] = (y - SY / 2) * 3.0f; q[2] = -30.0f;
    }
    for (int z = 0; z < DZ; ++z) for (int y = 0; y < DY; ++y) for (int x = 0; x < DX; ++x) {
        float* q = dp + 3 * ((z * DY + y) * DX + x);
        q[0] = x - DX / 2 + 0.5f; q[1] = y - DY / 2 + 0.5f; q[2] = -35.0f - z;
    }
    md_render_batch_t bs = {0}, bd = {0};
    bs.key = 1; bs.kind = MD_RENDER_SPHERES; bs.count = NS; bs.points.count = NS;
    bs.points.position = f_buf(&fx, sp, NS * 12);
    bs.style.spheres.radius_add = 0.5f;
    bd = bs;
    bd.count = ND; bd.points.count = ND;
    bd.points.position = f_buf(&fx, dp, (size_t)ND * 12);
    bd.style.spheres.radius_add = 0.8f;
    md_render_frame_t fsp = frame_of(W, H, make_camera(&c), &bs, 1), fde = frame_of(W, H, make_camera(&c), &bd, 1);
    md_render_frame_t rsp = fsp, rde = fde;
    md_render_view_t ref_s = md_render_view_create(fx.r, "sparse, no occlusion"), ref_d = md_render_view_create(fx.r, "dense, no occlusion");
    md_render_debug_view(ref_s, MD_RENDER_DEBUG_NO_COMPUTE_RASTER | MD_RENDER_DEBUG_NO_OCCLUSION, 0.0f);
    md_render_debug_view(ref_d, MD_RENDER_DEBUG_NO_COMPUTE_RASTER | MD_RENDER_DEBUG_NO_OCCLUSION, 0.0f);
    ASSERT_TRUE(md_render_view_draw(fx.gfx, ref_s, &rsp, NULL));
    ASSERT_TRUE(md_render_view_draw(fx.gfx, ref_d, &rde, NULL));

    int hiz = 0, covered = 0;
    md_render_debug_pass_hook(count_hiz, &hiz);
    md_render_view_t v = md_render_view_create(fx.r, "occlusion as it pays");
    md_render_debug_view(v, MD_RENDER_DEBUG_NO_COMPUTE_RASTER, 0.0f);
    /* sparse: the first raster builds the pyramid, the second tests it and
       culls nothing, the ones after go without */
    EXPECT_EQ(1, raster_with_hiz(&fx, v, &fsp, &hiz));
    EXPECT_EQ(1, raster_with_hiz(&fx, v, &fsp, &hiz));
    EXPECT_EQ(0, raster_with_hiz(&fx, v, &fsp, &hiz));
    EXPECT_EQ(0, raster_with_hiz(&fx, v, &fsp, &hiz));
    EXPECT_EQ(0, compare_views(&fx, v, ref_s, W, H, &covered));
    EXPECT_GT(covered, 500);
    /* dense, in the same view: drawn without (it was off), then the elements
       drawn more than doubled, so on again: a pyramid, then culling with
       it, which pays, so it stays */
    EXPECT_EQ(0, raster_with_hiz(&fx, v, &fde, &hiz));
    EXPECT_EQ(1, raster_with_hiz(&fx, v, &fde, &hiz));
    EXPECT_EQ(1, raster_with_hiz(&fx, v, &fde, &hiz));
    md_render_stats_t st;
    ASSERT_TRUE(md_render_view_stats(v, &st));
    EXPECT_GT(st.occluded, (uint32_t)ND / 2);
    EXPECT_EQ(1, raster_with_hiz(&fx, v, &fde, &hiz));
    EXPECT_EQ(1, raster_with_hiz(&fx, v, &fde, &hiz));
    EXPECT_EQ(0, compare_views(&fx, v, ref_d, W, H, &covered));
    EXPECT_GT(covered, (int)(W * H / 4));
    printf("    dense: %u of %u culled by occlusion, %u drawn\n", st.occluded, (unsigned)ND, st.hardware[0] + st.hardware[1]);
    /* the switch off: every raster culls, however little */
    md_render_debug_occlusion_auto(false);
    md_render_view_t va = md_render_view_create(fx.r, "occlusion always");
    md_render_debug_view(va, MD_RENDER_DEBUG_NO_COMPUTE_RASTER, 0.0f);
    for (int k = 0; k < 4; ++k) EXPECT_EQ(1, raster_with_hiz(&fx, va, &fsp, &hiz));
    md_render_debug_occlusion_auto(true);
    md_render_debug_pass_hook(NULL, NULL);
    md_render_view_destroy(fx.gfx, va);
    md_render_view_destroy(fx.gfx, v);
    md_render_view_destroy(fx.gfx, ref_s);
    md_render_view_destroy(fx.gfx, ref_d);
    free(sp); free(dp);
    f_close(&fx);
}

/* Radii from tables: six tables over seven batches (one read by two, one of
   600 entries), types past a table's end among them, small spheres for the compute raster and large ones for the hardware. They
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
            memcpy(&pts[4 * i + 3], &type, 4);
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
    if (info.atomic_int64_min_max) EXPECT_GT(st.compute[0] + st.compute[1], 0u);
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
    md_gpu_addr_t P = f_buf(&fx, pos, N * 12), R = f_buf(&fx, rad, N * 4), S = f_marks(&fx, sel, NULL, N);
    fs.vis[0] = f_buf(&fx, NULL, N / 8);
    fs.vis[1] = f_buf(&fx, NULL, N / 8);
    md_render_camera_t cam = make_camera(&c);
    const float px = fabsf(cam.proj[5]) * 0.5f * (float)H;
    /* the band: from sphere 150's depth to sphere 300's, inside clusters */
    float band[2] = { px / -fs.z[150], px / -fs.z[300] };
    for (int k = 0; k < 2; ++k) {
        md_render_batch_t b = {0};
        b.key = 1 + k; b.kind = MD_RENDER_SPHERES; b.count = N; b.points.count = N;
        b.points.position = P; b.points.radius = R; b.points.marks = S;
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

/* Keyed THIN: groups of spheres, each with its centroid as the key of its
   members, fading out across a band, and the centroids themselves (a batch
   over the tail of the same points, its clusters the tail of the same
   offsets) fading in across it; sticks within the groups keyed alike. They
   draw what the same batches draw unfaded with each group given by
   visibility bits to the side its centroid's hash and depth put it on:
   every group as its members or as its stand-in, never both or neither, to
   the bit. A camera move hands groups over without a refit. */
typedef struct key_scene_t {
    uint32_t      groups, members;
    float*        cz;              /* view-space z of each centroid, camera at rest */
    md_gpu_addr_t vis[2];          /* the reference's bits over the attributes: members, stand-ins */
} key_scene_t;

static double key_split(fix_t* fx, key_scene_t* ks, const float band[2], float px, float dz, uint32_t* to_standin) {
    const uint32_t G = ks->groups, M = ks->members, A = G * M + G, words = (A + 31) / 32;
    uint32_t* bits = calloc(2 * words, 4);
    float ab[2];
    fade_ab(band, px, false, ab);
    double margin = 1e30;
    *to_standin = 0;
    for (uint32_t g = 0; g < G; ++g) {
        const uint32_t sp = G * M + g;              /* the stand-in's point and attribute */
        float d = -(ks->cz[g] + dz);
        float t = fade_t_cpu(ab, d), h = fade_h_cpu(sp);
        if (h >= t) { for (uint32_t m = 0; m < M; ++m) bits[(g * M + m) / 32] |= 1u << ((g * M + m) % 32); }
        else        { bits[words + sp / 32] |= 1u << (sp % 32); (*to_standin)++; }
        if (t > 0.0f && t < 1.0f && fabs((double)h - t) < margin) margin = fabs((double)h - t);
    }
    md_gpu_upload(fx->gfx, ks->vis[0], bits, words * 4);
    md_gpu_upload(fx->gfx, ks->vis[1], bits + words, words * 4);
    free(bits);
    return margin;
}

UTEST(render, fade_key_hands_over_groups) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 160, H = 120;
    cam_t c = { .tan_half = 0.6f, .aspect = (float)W / H, .n = 0.1f };
    enum { G = 96, M = 8, NP = G * M + G, NS = G * (M - 1), PER = 8 };   /* PER groups a cluster */
    key_scene_t ks = { .groups = G, .members = M };
    ks.cz = malloc(G * sizeof(float));
    float*    pos  = malloc(NP * 12);
    float*    rad  = malloc(NP * 4);
    uint32_t* key  = malloc(NP * 4);
    uint32_t* attr = malloc(G * 4);
    uint32_t* pair = malloc(NS * 8);
    uint32_t  coff[G / PER + G / PER + 1];          /* the members' clusters, then the stand-ins' */
    uint32_t  soff[G / PER + 1];                    /* the sticks' */
    uint32_t  sel[(NP + 31) / 32];
    memset(sel, 0x55, sizeof(sel));                 /* every other attribute selected: masks */
    uint32_t seed = 7;
    #define RND() ((seed = seed * 1664525u + 1013904223u) >> 8) / 16777216.0f
    for (uint32_t g = 0; g < G; ++g) {              /* farther with g */
        float z = -4.0f - 40.0f * (float)g / G;
        float cx = (RND() * 2.0f - 1.0f) * 0.5f * -z * c.tan_half * c.aspect;
        float cy = (RND() * 2.0f - 1.0f) * 0.5f * -z * c.tan_half;
        double sum[3] = { 0, 0, 0 };
        for (uint32_t m = 0; m < M; ++m) {
            uint32_t i = g * M + m;
            pos[3 * i + 0] = cx + 0.4f * (RND() - 0.5f);
            pos[3 * i + 1] = cy + 0.4f * (RND() - 0.5f);
            pos[3 * i + 2] = z  + 0.4f * (RND() - 0.5f);
            rad[i] = 0.06f + 0.06f * RND();
            for (int k = 0; k < 3; ++k) sum[k] += pos[3 * i + k];
            key[i] = G * M + g;
        }
        uint32_t sp = G * M + g;
        for (int k = 0; k < 3; ++k) pos[3 * sp + k] = (float)(sum[k] / M);
        rad[sp] = 0.25f;
        key[sp] = sp;
        attr[g] = sp;
        ks.cz[g] = pos[3 * sp + 2];
        for (uint32_t m = 0; m + 1 < M; ++m) { pair[2 * (g * (M - 1) + m) + 0] = g * M + m; pair[2 * (g * (M - 1) + m) + 1] = g * M + m + 1; }
    }
    #undef RND
    for (uint32_t k = 0; k <= G / PER; ++k) { coff[k] = k * PER * M; soff[k] = k * PER * (M - 1); }
    for (uint32_t k = 1; k <= G / PER; ++k) coff[G / PER + k] = G * M + k * PER;   /* coff[G / PER] = G M: the stand-ins' first */
    md_gpu_addr_t P = f_buf(&fx, pos, NP * 12), R = f_buf(&fx, rad, NP * 4), K = f_buf(&fx, key, NP * 4);
    md_gpu_addr_t AI = f_buf(&fx, attr, G * 4), PR = f_buf(&fx, pair, NS * 8), S = f_marks(&fx, sel, NULL, (uint32_t)sizeof(sel) * 8);
    md_gpu_addr_t CO = f_buf(&fx, coff, sizeof(coff)), SO = f_buf(&fx, soff, sizeof(soff));
    ks.vis[0] = f_buf(&fx, NULL, sizeof(sel));
    ks.vis[1] = f_buf(&fx, NULL, sizeof(sel));
    md_render_camera_t cam = make_camera(&c);
    const float px = fabsf(cam.proj[5]) * 0.5f * (float)H;
    float band[2] = { px / -ks.cz[20], px / -ks.cz[70] };
    uint32_t n_in = 0;
    while (key_split(&fx, &ks, band, px, 0.0f, &n_in) < 1e-4) band[1] *= 0.9993f;

    md_render_batch_t b[3], r[3];
    memset(b, 0, sizeof(b));
    for (int k = 0; k < 3; ++k) {
        b[k].key = 1 + k;
        b[k].points.position = P; b[k].points.radius = R; b[k].points.count = NP; b[k].points.attr_count = NP;
        b[k].points.marks = S;
        b[k].color.mode = MD_RENDER_COLOR_UNIFORM;
        b[k].color.uniform = 0xFF0000FFu + 0x4000u * k;
    }
    /* 0: the members, keyed, fading out */
    b[0].kind = MD_RENDER_SPHERES; b[0].count = G * M;
    b[0].cluster_offset = CO; b[0].cluster_count = G / PER;
    b[0].style.spheres.radius_scale = 1.0f;
    /* 1: the stand-ins, the tail of the points and of the cluster offsets */
    b[1].kind = MD_RENDER_SPHERES; b[1].count = G;
    b[1].points.position = P + G * M * 12; b[1].points.radius = R + G * M * 4; b[1].points.count = G;
    b[1].points.attr_index = AI;
    b[1].cluster_offset = CO + (G / PER) * 4; b[1].cluster_count = G / PER;
    b[1].style.spheres.radius_scale = 1.0f;
    /* 2: sticks within the groups, keyed */
    b[2].kind = MD_RENDER_CAPSULES; b[2].count = NS; b[2].index = PR;
    b[2].cluster_offset = SO; b[2].cluster_count = G / PER;
    b[2].style.capsules.radius = 0.03f;
    memcpy(r, b, sizeof(b));
    r[0].points.visible[0] = r[2].points.visible[0] = ks.vis[0];
    r[1].points.visible[0] = ks.vis[1];
    for (int k = 0; k < 3; ++k) {
        b[k].fade.mode = MD_RENDER_FADE_THIN;
        memcpy(k == 1 ? b[k].fade.in : b[k].fade.out, band, sizeof(band));
        if (k != 1) b[k].fade.key = K;
    }

    md_render_view_t ref = md_render_view_create(fx.r, "visibility bits");
    md_render_frame_t f  = frame_of(W, H, cam, b, 3);
    md_render_frame_t fr = frame_of(W, H, cam, r, 3);
    md_render_stats_t st;
    int covered = 0, set = 0;
    for (int pass = 0; pass < 2; ++pass) {
        if (pass == 1) {                              /* the camera moves in: groups change sides */
            float dz = 5.0f;
            while (key_split(&fx, &ks, band, px, dz, &n_in) < 1e-4) dz += 0.0137f;
            f.camera.view[14] = fr.camera.view[14] = dz;
            fr.content[MD_RENDER_CONTENT_GEOMETRY]++;     /* the reference's bits changed */
        }
        ASSERT_TRUE(md_render_view_draw(fx.gfx, ref, &fr, NULL));
        ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
        md_gpu_stream_sync(fx.gfx);
        ASSERT_TRUE(md_render_view_stats(fx.view, &st));
        if (pass == 1) EXPECT_FALSE(st.refitted);
        /* each group faded once: as members (and their sticks) or as its stand-in */
        EXPECT_EQ(n_in * (M + (M - 1)) + (G - n_in), st.faded);
        EXPECT_GT(n_in, 0u);
        EXPECT_LT(n_in, (uint32_t)G);
        EXPECT_EQ(0, compare_views(&fx, fx.view, ref, W, H, &covered));
        EXPECT_GT(covered, (int)(W * H / 50));
        EXPECT_EQ(0, compare_masks(&fx, fx.view, ref, W, H, &set));
        EXPECT_GT(set, 0);
        printf("    %u groups of %u, %u as their stand-in: %d pixels covered, %d masked, faded %u\n",
               (unsigned)G, (unsigned)M, n_in, covered, set, st.faded);
    }
    md_render_view_destroy(fx.gfx, ref);
    free(pos); free(rad); free(key); free(attr); free(pair); free(ks.cz);
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
    ASSERT_TRUE(md_system_infer(&sys, &state, NULL));
    ASSERT_EQ(8u, (uint32_t)(sys.bond.count + sys.bond.interaction_count));

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
    EXPECT_EQ(8u, topo.bond_kind_offset[1]);   /* all covalent */

    /* The caller's atom arrays */
    md_gpu_addr_t P = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, 9 * 12).gpu;
    ASSERT_TRUE(md_render_mol_upload_state(fx.gfx, P, &state, content));
    uint32_t colors[9];
    for (int i = 0; i < 9; ++i) colors[i] = 0xFF000000u | (uint32_t)(i * 20);
    md_render_mol_atoms_t atoms = {0};
    atoms.position = P;
    md_render_mol_ctx_t kernels = md_render_mol_create(fx.dev, NULL);
    ASSERT_TRUE(kernels != NULL);
    md_gpu_addr_t dyn = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, topo.dynamic_size).gpu;
    ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, 0, content));   /* the copy that is drawn */

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
    md_render_mol_destroy(kernels);
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
    const uint32_t nb_sys = (uint32_t)(sys->bond.count + sys->bond.interaction_count);
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
    for (uint32_t c = 0; c < MD_RENDER_MOL_BOND_KIND_COUNT; ++c) {
        const uint32_t* bo = (const uint32_t*)(img + T->bond_cluster_offset[c]);
        uint32_t size = T->bond_kind_offset[c + 1] - T->bond_kind_offset[c], cc = T->bond_cluster_count[c];
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
    return ok && md_system_infer(sys, state, NULL);
}

static uint8_t* read_block(fix_t* fx, md_gpu_addr_t base, uint64_t size);

/* The plan of systems with chains, ligands and small molecules; then (with a
   GPU) the same after a regroup, which moves the small molecules and with
   them the bonds that reach them from elsewhere (mixed) */
UTEST(render, mol_plan_clusters) {
    const char* files[] = { MD_UNITTEST_DATA_DIR "/tubulin-A-B.pdb", MD_UNITTEST_DATA_DIR "/1a64.pdb", MD_UNITTEST_DATA_DIR "/tryptophan-md.gro" };
    fix_t fx;
    bool gpu = f_open(&fx);
    md_render_mol_ctx_t kernels = gpu ? md_render_mol_create(fx.dev, NULL) : NULL;
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
        if (kernels) {
            md_render_mol_topo_t topo;
            ASSERT_TRUE(md_render_mol_write(fx.gfx, &topo, &plan, md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, plan.size).gpu, NULL));
            md_gpu_addr_t P   = f_buf(&fx, state.xyz, T->atom_count * 12);
            md_gpu_addr_t dyn = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, topo.dynamic_size).gpu;
            md_render_mol_atoms_t atoms = {0};
            atoms.position = P;
            md_render_mol_topo_t rel = topo;
            rel.atom_order -= topo.base; rel.mixed_point -= topo.base; rel.atom_cluster_offset -= topo.base;
            rel.bond_pairs -= topo.base; rel.bond_ids -= topo.base;
            for (int c = 0; c < MD_RENDER_MOL_BOND_KIND_COUNT; ++c) rel.bond_cluster_offset[c] -= topo.base;
            /* twice: the second from the first's records (the other buffer),
               mixed bonds and all */
            for (int rg = 0; rg < 2; ++rg) {
                ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, MD_RENDER_MOL_UPDATE_REGROUP, NULL));
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
    if (kernels) md_render_mol_destroy(kernels);
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

/* Spacefill's level of detail. The residues' points of a protein with its
   crystal waters:
   each run of a residue's static points in one cluster, the key of each of
   its atoms; their centroids and radii (sqrt(5/3 Rg^2 + mean r^2)) as the
   update writes them. Emitted with a band: the static atoms keyed, the
   water's atoms apart, the residues' spheres; drawn, they are the same
   batches unfaded with each residue given by visibility bits to its atoms
   or to its sphere by its point's depth and hash, to the bit, every residue
   once. A pick on a sphere is its residue's. Seen from far, the spheres
   cover about what their atoms do. */
static int count_covered(fix_t* fx, md_render_view_t v, uint32_t W, uint32_t H) {
    uint32_t* id = malloc(W * H * 8);
    int n = -1;
    if (f_read(fx, md_render_view_targets(v).id, id, W * H * 8)) { n = 0; for (uint32_t i = 0; i < W * H; ++i) n += id[2 * i + 1] != 0; }
    free(id);
    return n;
}

UTEST(render, mol_spacefill_residues) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(load_system(&sys, &state, MD_UNITTEST_DATA_DIR "/1a64.pdb"));
    const uint32_t n = (uint32_t)sys.atom.count;
    uint64_t content[MD_RENDER_CONTENT_COUNT] = {0};

    md_temp_scope_t temp = md_temp_begin();
    md_render_mol_plan_t plan;
    ASSERT_TRUE(md_render_mol_plan(&plan, &sys, NULL, temp));
    const md_render_mol_topo_t* T = &plan.layout;
    const uint32_t n_static = T->point_count - T->small_atom_count, n_rp = T->residue_point_count, first = T->point_count + T->derived_count;
    ASSERT_GT(n_rp, 0u);
    ASSERT_GT(T->small_atom_count, 0u);
    uint32_t* order = malloc((first + n_rp) * 4);   /* host copies: the image goes with the scope */
    uint32_t* rp_off = malloc((n_rp + 1) * 4), *rp_res = malloc(n_rp * 4), *key = malloc(n_static * 4);
    memcpy(order,  (const uint8_t*)plan.data + T->atom_order, (first + n_rp) * 4);
    memcpy(rp_off, (const uint8_t*)plan.data + T->residue_point_offset, (n_rp + 1) * 4);
    memcpy(rp_res, (const uint8_t*)plan.data + T->residue_point_residue, n_rp * 4);
    memcpy(key,    (const uint8_t*)plan.data + T->residue_key, n_static * 4);
    {
        int bad = 0;
        const uint32_t* coff = (const uint32_t*)((const uint8_t*)plan.data + T->atom_cluster_offset);
        const uint32_t* rcof = (const uint32_t*)((const uint8_t*)plan.data + T->residue_cluster_offset);
        uint32_t c = 0;
        EXPECT_EQ(0u, rp_off[0]);
        EXPECT_EQ(n_static, rp_off[n_rp]);
        for (uint32_t k = 0; k < n_rp && bad < 5; ++k) {
            const uint32_t p0 = rp_off[k], p1 = rp_off[k + 1], ra = order[first + k];
            if (p1 <= p0) { bad++; printf("    residue point %u empty\n", k); continue; }
            while (coff[c + 1] <= p0) ++c;
            if (p1 > coff[c + 1]) { bad++; printf("    residue point %u crosses cluster %u\n", k, c); }
            bool has = false;
            for (uint32_t q = p0; q < p1; ++q) {
                const uint32_t a = order[q];
                const uint32_t r = rp_res[k];
                if (r != UINT32_MAX && !(a >= sys.residue.atom_offset[r] && a < sys.residue.atom_offset[r + 1])) { bad++; printf("    atom %u not of residue %u\n", a, r); }
                if (key[q] != first + k) { bad++; printf("    point %u keyed %u, not %u\n", q, key[q], first + k); }
                has |= a == ra;
            }
            if (!has) { bad++; printf("    residue point %u stands for atom %u, not its own\n", k, ra); }
        }
        EXPECT_EQ(0, bad);
        EXPECT_EQ(0u, rcof[0]);
        EXPECT_EQ(n_rp, rcof[T->residue_cluster_count]);
        for (uint32_t k = 0; k < T->residue_cluster_count; ++k) EXPECT_LE(rcof[k + 1] - rcof[k], 64u);
        printf("    %u atoms, %u static: %u residue points in %u clusters\n", n, n_static, n_rp, T->residue_cluster_count);
    }
    md_gpu_mem_t mem = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, plan.size);
    md_render_mol_topo_t topo;
    ASSERT_TRUE(md_render_mol_write(fx.gfx, &topo, &plan, mem.gpu, content));
    md_temp_end(temp);

    md_gpu_addr_t P = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, (size_t)n * 12).gpu;
    ASSERT_TRUE(md_render_mol_upload_state(fx.gfx, P, &state, content));
    md_render_mol_atoms_t atoms = {0};
    atoms.position = P;
    md_render_mol_ctx_t kernels = md_render_mol_create(fx.dev, NULL);
    ASSERT_TRUE(kernels != NULL);
    md_gpu_addr_t dyn = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, topo.dynamic_size).gpu;
    ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, 0, content));
    md_gpu_stream_sync(fx.gfx);

    /* The update's points against the host's */
    float* rp = malloc(n_rp * 16);
    {
        uint8_t* blk = read_block(&fx, dyn + (uint64_t)first * 16u, (uint64_t)n_rp * 16u);
        ASSERT_TRUE(blk != NULL);
        memcpy(rp, blk, n_rp * 16);
        free(blk);
        double worst_c = 0, worst_r = 0;
        for (uint32_t k = 0; k < n_rp; ++k) {
            const uint32_t p0 = rp_off[k], p1 = rp_off[k + 1], m = p1 - p0;
            double s[3] = { 0, 0, 0 }, s2 = 0, r2 = 0;
            const float* ref = &state.xyz[order[p0]].x;
            for (uint32_t q = p0; q < p1; ++q) {
                const float* x = &state.xyz[order[q]].x;
                double d[3] = { x[0] - ref[0], x[1] - ref[1], x[2] - ref[2] };
                for (int j = 0; j < 3; ++j) s[j] += d[j];
                s2 += d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                const double r = sys.atom.type.radius[sys.atom.type_idx[order[q]]];
                r2 += r * r;
            }
            double c[3], cc = 0;
            for (int j = 0; j < 3; ++j) { c[j] = s[j] / m; cc += c[j] * c[j]; }
            double rad = sqrt(5.0 / 3.0 * fmax(s2 / m - cc, 0.0) + r2 / m);
            for (int j = 0; j < 3; ++j) worst_c = fmax(worst_c, fabs(ref[j] + c[j] - rp[4 * k + j]));
            worst_r = fmax(worst_r, fabs(rad - rp[4 * k + 3]));
        }
        EXPECT_LT(worst_c, 1e-4);
        EXPECT_LT(worst_r, 1e-4);
        printf("    centroids within %.2g A, radii within %.2g A of the host's\n", worst_c, worst_r);
    }

    md_render_mol_rep_t rep = {0};
    rep.key  = 7;
    rep.type = MD_RENDER_MOL_REP_SPACEFILL;
    rep.params.spacefill.radius_scale = 1.0f;
    rep.color.mode = MD_RENDER_COLOR_UNIFORM;
    rep.color.uniform = 0xFF8080FFu;
    md_render_batch_t b[4];
    ASSERT_EQ(1u, md_render_mol_emit(&topo, dyn, &atoms, &rep, b, 4));      /* no band: the atoms */
    md_render_batch_t atoms_only = b[0];

    /* A perspective view down -z at the protein from outside the box */
    const uint32_t W = 160, H = 120;
    cam_t c = { .tan_half = 0.5f, .aspect = (float)W / H, .n = 1.0f };
    md_render_camera_t cam = make_camera(&c);
    double ctr[3] = { 0, 0, 0 }, zmax = -1e30;
    for (uint32_t q = 0; q < n_static; ++q) for (int j = 0; j < 3; ++j) ctr[j] += (&state.xyz[order[q]].x)[j] / n_static;
    for (uint32_t i = 0; i < n; ++i) zmax = fmax(zmax, state.xyz[i].z);
    cam.view[12] = (float)-ctr[0]; cam.view[13] = (float)-ctr[1]; cam.view[14] = (float)-(zmax + 40.0);
    const float px = fabsf(cam.proj[5]) * 0.5f * (float)H;
    /* the residues' view depths; the band over their middle */
    float* d = malloc(n_rp * 4);
    float dmin = 1e30f, dmax = 0;
    for (uint32_t k = 0; k < n_rp; ++k) {
        d[k] = -(rp[4 * k + 2] + cam.view[14]);
        dmin = fminf(dmin, d[k]); dmax = fmaxf(dmax, d[k]);
    }
    float band[2] = { px / (dmin + 0.3f * (dmax - dmin)), px / (dmin + 0.7f * (dmax - dmin)) };
    const uint32_t words = (n + 31) / 32;
    uint32_t* bits = malloc(2 * words * 4);
    md_gpu_addr_t VA = f_buf(&fx, NULL, words * 4), VR = f_buf(&fx, NULL, words * 4);
    uint32_t to_res = 0, faded = 0;
    for (int tries = 0; ; ++tries) {
        float ab[2];
        fade_ab(band, px, false, ab);
        double margin = 1e30;
        memset(bits, 0, 2 * words * 4);
        for (uint32_t i = 0; i < n; ++i) bits[i / 32] |= 1u << (i % 32);    /* the water's atoms, and the static ones kept */
        to_res = faded = 0;
        for (uint32_t k = 0; k < n_rp; ++k) {
            const uint32_t ra = order[first + k];
            const float t = fade_t_cpu(ab, d[k]), h = fade_h_cpu(ra);
            if (t > 0.0f && t < 1.0f) margin = fmin(margin, fabs((double)h - t));
            if (h < t) {
                to_res++;
                faded += rp_off[k + 1] - rp_off[k];
                for (uint32_t q = rp_off[k]; q < rp_off[k + 1]; ++q) bits[order[q] / 32] &= ~(1u << (order[q] % 32));
                bits[words + ra / 32] |= 1u << (ra % 32);
            } else {
                faded++;
            }
        }
        if (margin >= 1e-4 || tries > 200) break;
        band[1] *= 0.9993f;
    }
    md_gpu_upload(fx.gfx, VA, bits, words * 4);
    md_gpu_upload(fx.gfx, VR, bits + words, words * 4);
    memcpy(rep.params.spacefill.residues, band, sizeof(band));
    const uint32_t nb = md_render_mol_emit(&topo, dyn, &atoms, &rep, b, 4);
    ASSERT_EQ(3u, nb);
    EXPECT_EQ((uint64_t)7 << 8 | MD_RENDER_MOL_PART_ATOMS, b[0].key);
    EXPECT_EQ((uint64_t)7 << 8 | MD_RENDER_MOL_PART_SMALL_ATOMS, b[1].key);
    EXPECT_EQ((uint64_t)7 << 8 | MD_RENDER_MOL_PART_RESIDUES, b[2].key);
    EXPECT_EQ(n_static, b[0].count);
    EXPECT_EQ(T->small_atom_count, b[1].count);
    EXPECT_EQ(n_rp, b[2].count);
    EXPECT_EQ(topo.residue_key, b[0].fade.key);
    EXPECT_EQ(0u, b[1].fade.mode);
    md_render_batch_t r[3];
    memcpy(r, b, sizeof(r));
    for (int k = 0; k < 3; ++k) memset(&r[k].fade, 0, sizeof(r[k].fade));
    r[0].points.visible[0] = r[1].points.visible[0] = VA;
    r[2].points.visible[0] = VR;

    md_render_view_t ref = md_render_view_create(fx.r, "visibility bits");
    md_render_frame_t f  = frame_of(W, H, cam, b, nb);
    md_render_frame_t fr = frame_of(W, H, cam, r, nb);
    memcpy(f.content, content, sizeof(content));
    memcpy(fr.content, content, sizeof(content));
    ASSERT_TRUE(md_render_view_draw(fx.gfx, ref, &fr, NULL));
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    md_gpu_stream_sync(fx.gfx);
    md_render_stats_t st;
    ASSERT_TRUE(md_render_view_stats(fx.view, &st));
    int covered = 0;
    EXPECT_EQ(0, compare_views(&fx, fx.view, ref, W, H, &covered));
    EXPECT_EQ(faded, st.faded);
    EXPECT_GT(to_res, n_rp / 10);
    EXPECT_LT(to_res, n_rp - n_rp / 10);
    printf("    band %.2f .. %.2f px/A: %u of %u residues as their sphere, %d pixels covered, faded %u\n",
           band[0], band[1], to_res, n_rp, covered, st.faded);

    /* All residues as spheres, the water hidden: a pick is a residue's, and
       from far the spheres cover about what the atoms do */
    rep.params.spacefill.residues[0] = 1e6f; rep.params.spacefill.residues[1] = 1e5f;
    ASSERT_EQ(3u, md_render_mol_emit(&topo, dyn, &atoms, &rep, b, 4));
    md_render_batch_t only[2] = { b[2], atoms_only };
    only[1].count = n_static;
    only[1].cluster_count = topo.atom_cluster_count - topo.grid_cells;
    cam.view[14] = (float)-(zmax + 250.0);
    for (int k = 0; k < 2; ++k) {
        md_render_frame_t fo = frame_of(W, H, cam, &only[k], 1);
        memcpy(fo.content, content, sizeof(content));
        fo.content[MD_RENDER_CONTENT_GEOMETRY] += 1 + k;
        ASSERT_TRUE(md_render_view_draw(fx.gfx, k ? ref : fx.view, &fo, NULL));
    }
    md_gpu_stream_sync(fx.gfx);
    const int cov_res = count_covered(&fx, fx.view, W, H), cov_atoms = count_covered(&fx, ref, W, H);
    EXPECT_GT(cov_atoms, 200);
    EXPECT_NEAR(1.0, (double)cov_res / cov_atoms, 0.15);
    printf("    from %.0f A: the residues' spheres cover %d pixels, their atoms %d\n", zmax + 250.0 - ctr[2], cov_res, cov_atoms);
    {
        uint32_t* id = malloc(W * H * 8);
        ASSERT_TRUE(f_read(&fx, md_render_view_targets(fx.view).id, id, W * H * 8));
        uint32_t picked = 0, seen = 0;
        for (uint32_t i = 0; i < W * H && picked < 8; ++i) {
            if (!id[2 * i + 1] || seen++ % 31) continue;           /* every 31st covered pixel */
            md_render_pick_t pk;
            md_render_mol_hit_t hit;
            ASSERT_TRUE(f_pick(&fx, i % W, i / W, &pk));
            ASSERT_TRUE(md_render_mol_decode_pick(&pk, &hit));
            EXPECT_EQ(MD_RENDER_MOL_HIT_RESIDUE, hit.kind);
            ASSERT_LT(hit.residue, (uint32_t)sys.residue.count);
            EXPECT_GE(hit.atom, sys.residue.atom_offset[hit.residue]);
            EXPECT_LT(hit.atom, sys.residue.atom_offset[hit.residue + 1]);
            picked++;
        }
        free(id);
        EXPECT_EQ(8u, picked);
    }

    md_render_view_destroy(fx.gfx, ref);
    free(order); free(rp_off); free(rp_res); free(key); free(rp); free(d); free(bits);
    md_gpu_free(fx.gfx, P);
    md_gpu_free(fx.gfx, dyn);
    md_gpu_free(fx.gfx, topo.base);
    md_render_mol_destroy(kernels);
    md_system_free(&sys);
    f_close(&fx);
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
    md_render_mol_ctx_t kernels = md_render_mol_create(fx.dev, NULL);
    ASSERT_TRUE(kernels != NULL);
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
        ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, steps[k].atoms, steps[k].flags, NULL));
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
        ASSERT_EQ(pa.layout.point_count, pb.layout.point_count);   /* (the sizes may differ: the bonds left out differ) */
        md_gpu_addr_t base = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, MD_MAX(pa.size, pb.size)).gpu;
        md_gpu_addr_t dyn2 = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, MD_MAX(pa.layout.dynamic_size, pb.layout.dynamic_size)).gpu;
        md_render_mol_topo_t ta, tb;
        ASSERT_TRUE(md_render_mol_write(fx.gfx, &ta, &pa, base, NULL));
        ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &ta, dyn2, &types, 0, NULL));
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
            ASSERT_TRUE(md_render_mol_frame_upload(kernels, fx.gfx, &ta, dyn2, &fr, 0, &f, 0));
            md_render_mol_atoms_t af = {0};
            af.frames = &fr; af.frame_count = 1;
            EXPECT_TRUE(md_render_mol_update(kernels, fx.gfx, &ta, dyn2, &af, 0, NULL));
            md_gpu_stream_sync(fx.gfx);
            free(q);
        }
        ASSERT_TRUE(md_render_mol_write(fx.gfx, &tb, &pb, base, NULL));
        EXPECT_NE(ta.instance, tb.instance);
        ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &tb, dyn2, &types, 0, NULL));
        const int wrong = radii_wrong(&fx, &tb, dyn2, &types, &rep, by_type);
        if (wrong) printf("    written again at the same base: %d points with a wrong radius\n", wrong);
        EXPECT_EQ(0, wrong);
        md_render_mol_atoms_t af = {0};
        af.frames = &fr; af.frame_count = 1;
        EXPECT_FALSE(md_render_mol_update(kernels, fx.gfx, &tb, dyn2, &af, 0, NULL));   /* frames of the other one */
        md_temp_end(t2);
        md_gpu_free(fx.gfx, frb); md_gpu_free(fx.gfx, dyn2); md_gpu_free(fx.gfx, base);
        md_bitmask_free(&ra); md_bitmask_free(&rb);
    }
    free(by_type); free(own);
    md_gpu_free(fx.gfx, types.position); md_gpu_free(fx.gfx, mine.radius);
    md_gpu_free(fx.gfx, dyn); md_gpu_free(fx.gfx, topo.base);
    md_render_mol_destroy(kernels);
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
    const uint32_t n = (uint32_t)sys.atom.count, nb = (uint32_t)(sys.bond.count + sys.bond.interaction_count);
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

    md_render_mol_ctx_t kernels = md_render_mol_create(fx.dev, NULL);
    ASSERT_TRUE(kernels != NULL);
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
    ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, MD_RENDER_MOL_UPDATE_REGROUP, content));
    uint8_t* after = read_block(&fx, topo.base, topo.size);
    /* the block's own image: rebase the addresses to offsets */
    md_render_mol_topo_t rel = topo;
    rel.atom_order -= topo.base; rel.mixed_point -= topo.base; rel.atom_cluster_offset -= topo.base; rel.bond_pairs -= topo.base; rel.bond_ids -= topo.base;
    for (int c = 0; c < MD_RENDER_MOL_BOND_KIND_COUNT; ++c) rel.bond_cluster_offset[c] -= topo.base;
    EXPECT_EQ(0, check_topo_image(&sys, &rel, after, NULL));
    double e0 = grid_cluster_extent(&rel, before, xyz), e1 = grid_cluster_extent(&rel, after, xyz);
    printf("    grid clusters: mean extent %.1f A in molecule order, %.1f A after the update\n", e0, e1);
    EXPECT_LT(e1, 0.4 * e0);
    /* again from the regrouped order: the same cells (the order within one is the atomics') */
    ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, MD_RENDER_MOL_UPDATE_REGROUP, content));
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
        ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn_v, &av, 0, content));
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
    rep.bond_kinds = MD_RENDER_MOL_BONDS_ALL;
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
    md_render_debug_view(ref, MD_RENDER_DEBUG_NO_OCCLUSION | MD_RENDER_DEBUG_NO_COMPUTE_RASTER, 0.0f);
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
            ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, MD_RENDER_MOL_UPDATE_RADII, content));   /* refills the copy too */
        }
        /* a plain update (positions only) keeps the copy's radii: its w, not the source */
        ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, 0, content));
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
    md_render_mol_destroy(kernels);
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
    md_render_mol_ctx_t kernels = md_render_mol_create(fx.dev, NULL);
    ASSERT_TRUE(kernels != NULL);
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
    ASSERT_TRUE(md_render_mol_frame_upload(kernels, fx.gfx, &topo, dyn, &fr, 0, &frames[0], RG));
    uint8_t* order1 = read_block(&fx, topo.atom_order, (size_t)nres * 4);
    ASSERT_TRUE(md_render_mol_frame_upload(kernels, fx.gfx, &topo, 0,   &fr, 1, &frames[1], 0));
    ASSERT_TRUE(md_render_mol_frame_upload(kernels, fx.gfx, &topo, 0,   &fr, 2, &frames[2], 0));
    ASSERT_TRUE(md_render_mol_frame_upload(kernels, fx.gfx, &topo, dyn, &fr, 3, &frames[3], RG));
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
            ASSERT_TRUE(md_render_mol_frame_upload(kernels, fx.gfx, &topo, dyn, &fr, 0, &frames[4], RG));
            held[0] = 4;
        }
        atoms.frame_count = cases[ci].count;
        memcpy(atoms.frame_slot, cases[ci].slot, sizeof(atoms.frame_slot));
        md_render_mol_frame_weights(atoms.frame_w, atoms.frame_dw, cases[ci].count, cases[ci].t, cases[ci].dt, 0.5f);
        ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, 0, NULL));
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
            ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, 0, NULL));
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
    for (int c = 0; c < MD_RENDER_MOL_BOND_KIND_COUNT; ++c) rel.bond_cluster_offset[c] -= topo.base;
    EXPECT_EQ(0, check_topo_image(&sys, &rel, img, resident));
    /* refused: a regroup in the update, a slot without a frame, a cell with a
       not along x (the minimum image is GROMACS's, for that form) */
    EXPECT_FALSE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, RG, NULL));
    {
        md_render_mol_frame_t tilted = frames[0];
        tilted.cell[1] = 0.5f;
        EXPECT_FALSE(md_render_mol_frame_upload(kernels, fx.gfx, &topo, 0, &fr, 0, &tilted, 0));
    }
    md_render_mol_frames_t empty;
    ASSERT_TRUE(md_render_mol_frames_init(&empty, &topo, fr.base, 4));
    atoms.frames = &empty;
    EXPECT_FALSE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, 0, NULL));
    #undef RND

    md_gpu_stream_sync(fx.gfx);
    free(img); free(order0); free(order1); free(order2); free(v); free(xyz);
    for (int f = 0; f < NF; ++f) { free(fpos[f]); free(deq[f]); free(qs[f]); }
    md_gpu_free(fx.gfx, up[0].gpu);
    md_gpu_free(fx.gfx, up[1].gpu);
    md_gpu_free(fx.gfx, fr.base);
    md_gpu_free(fx.gfx, dyn);
    md_gpu_free(fx.gfx, topo.base);
    md_render_mol_destroy(kernels);
    md_system_free(&sys);
    md_bitmask_free(&mask);
    free(rlist);
    f_close(&fx);
}

/* Velocities derived from the positions (velocity_dt): each atom's
   displacement since the last update of the per-frame block, at the minimum
   image, over the time between, in the point order of the copy, through a
   regroup that moves the molecules a jump across the cell carried to other
   cells; at rest at a block's first update */
UTEST(render, mol_velocity_derived) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(load_system(&sys, &state, MD_UNITTEST_DATA_DIR "/tryptophan-md.gro"));
    const uint32_t n = (uint32_t)sys.atom.count;
    uint64_t content[MD_RENDER_CONTENT_COUNT] = {0};
    md_temp_scope_t temp = md_temp_begin();
    md_render_mol_plan_t plan;
    ASSERT_TRUE(md_render_mol_plan(&plan, &sys, NULL, temp));
    md_render_mol_topo_t topo;
    ASSERT_TRUE(md_render_mol_write(fx.gfx, &topo, &plan, md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, plan.size).gpu, content));
    md_temp_end(temp);
    ASSERT_GT(topo.small_molecule_count, 100u);
    md_render_mol_ctx_t kernels = md_render_mol_create(fx.dev, NULL);
    ASSERT_TRUE(kernels != NULL);
    md_gpu_addr_t dyn = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, topo.dynamic_size_velocity).gpu;

    float cell[9];
    md_render_mol_frame_cell(cell, &state.unitcell);
    ASSERT_GT(cell[0], 1.0f);
    float* x0 = malloc((size_t)n * 12);
    memcpy(x0, state.xyz, (size_t)n * 12);
    md_gpu_addr_t P = f_buf(&fx, x0, (size_t)n * 12);
    md_render_mol_atoms_t atoms = {0};
    atoms.position = P;
    memcpy(atoms.cell, cell, sizeof(cell));
    atoms.velocity_dt = 0.5f;
    md_render_mol_rep_t rep = {0};
    rep.key = 3; rep.type = MD_RENDER_MOL_REP_SPACEFILL; rep.params.spacefill.radius_scale = 1.0f;

    /* the first update of the block: at rest */
    ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, MD_RENDER_MOL_UPDATE_REGROUP, content));
    md_render_batch_t b;
    ASSERT_EQ(1u, md_render_mol_emit(&topo, dyn, &atoms, &rep, &b, 1));
    ASSERT_TRUE(b.points.velocity != 0);
    {
        uint8_t* blk = read_block(&fx, dyn, topo.dynamic_size_velocity);
        const uint32_t* v = (const uint32_t*)(blk + (b.points.velocity - dyn));
        int moving = 0;
        for (uint32_t p = 0; p < n; ++p) moving += (v[2 * p] & 0x7FFF7FFFu) != 0 || (v[2 * p + 1] & 0x7FFFu) != 0;
        EXPECT_EQ(0, moving);
        free(blk);
    }

    /* every atom moved by its own small displacement, every fifth small
       molecule as well across the cell (a jump by a and by -b) */
    float* d  = malloc((size_t)n * 12);
    float* x1 = malloc((size_t)n * 12);
    for (uint32_t i = 0; i < n; ++i) {
        d[3 * i + 0] =  0.20f * (float)(i % 7) / 7.0f;
        d[3 * i + 1] = -0.10f * (float)(i % 5) / 5.0f;
        d[3 * i + 2] =  0.05f * (float)(i % 3) / 3.0f;
        for (int k = 0; k < 3; ++k) x1[3 * i + k] = x0[3 * i + k] + d[3 * i + k];
    }
    const md_molecule_data_t* M = &sys.molecule;
    uint32_t jumped = 0;
    for (size_t m = 0; m < M->count; ++m) {
        const uint32_t sz = M->offset[m + 1] - M->offset[m];
        if (sz > MD_RENDER_MOL_SMALL || m % 5) continue;
        for (uint32_t s = M->offset[m]; s < M->offset[m + 1]; ++s) {
            float* x = x1 + 3 * M->atom_idx[s];
            for (int k = 0; k < 3; ++k) x[k] += cell[k] - cell[3 + k];
        }
        jumped++;
    }
    ASSERT_GT(jumped, 20u);
    md_gpu_upload(fx.gfx, P, x1, (size_t)n * 12);
    ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, MD_RENDER_MOL_UPDATE_REGROUP, content));
    ASSERT_EQ(1u, md_render_mol_emit(&topo, dyn, &atoms, &rep, &b, 1));
    {
        uint8_t* blk  = read_block(&fx, dyn, topo.dynamic_size_velocity);
        uint8_t* topb = read_block(&fx, topo.base, topo.size);
        const uint32_t* v     = (const uint32_t*)(blk + (b.points.velocity - dyn));
        const float*    pts   = (const float*)(blk + (b.points.position - dyn));
        const uint32_t* order = (const uint32_t*)(topb + (topo.atom_order - topo.base));
        int wrong = 0, misplaced = 0;
        double worst = 0.0;
        for (uint32_t p = 0; p < topo.point_count; ++p) {
            const uint32_t a = order[p];
            const float got[3] = { half_bits_to_float(v[2 * p] & 0xFFFFu), half_bits_to_float(v[2 * p] >> 16), half_bits_to_float(v[2 * p + 1] & 0xFFFFu) };
            for (int k = 0; k < 3; ++k) {
                const double want = d[3 * a + k] / 0.5, err = fabs(got[k] - want);
                worst = err > worst ? err : worst;
                wrong += err > 2e-3 + 1e-3 * fabs(want);
                misplaced += pts[4 * p + k] != x1[3 * a + k];
            }
        }
        printf("    velocities derived: %u molecules jumped, worst error %.2g\n", jumped, worst);
        EXPECT_EQ(0, wrong);
        EXPECT_EQ(0, misplaced);
        free(blk); free(topb);
    }

    /* a jump: at rest */
    ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, MD_RENDER_MOL_UPDATE_RESET_VELOCITY, content));
    {
        uint8_t* blk = read_block(&fx, dyn, topo.dynamic_size_velocity);
        const uint32_t* v = (const uint32_t*)(blk + (b.points.velocity - dyn));
        int moving = 0;
        for (uint32_t p = 0; p < n; ++p) moving += (v[2 * p] & 0x7FFF7FFFu) != 0 || (v[2 * p + 1] & 0x7FFFu) != 0;
        EXPECT_EQ(0, moving);
        free(blk);
    }

    /* one source of velocities at a time */
    md_render_mol_atoms_t both = atoms;
    both.velocity = P;
    EXPECT_FALSE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &both, 0, content));

    free(d); free(x1); free(x0);
    md_gpu_free(fx.gfx, P);
    md_gpu_free(fx.gfx, dyn);
    md_gpu_free(fx.gfx, topo.base);
    md_render_mol_destroy(kernels);
    md_vm_arena_destroy(alloc);
    f_close(&fx);
}

/* The caller's pairs (PAIRS) draw what the topology's bonds draw (BONDS)
   when they are those bonds: through the per-frame block's points (the
   index map from atoms), so through a regroup that moves molecules a jump
   carried to other cells, with the same derived velocities; a pick reports
   the pair's id */
UTEST(render, mol_pairs_draw_what_bonds_draw) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(load_system(&sys, &state, MD_UNITTEST_DATA_DIR "/tryptophan-md.gro"));
    const uint32_t n = (uint32_t)sys.atom.count;
    uint64_t content[MD_RENDER_CONTENT_COUNT] = {0};
    md_temp_scope_t temp = md_temp_begin();
    md_render_mol_plan_t plan;
    ASSERT_TRUE(md_render_mol_plan(&plan, &sys, NULL, temp));
    md_render_mol_topo_t topo;
    ASSERT_TRUE(md_render_mol_write(fx.gfx, &topo, &plan, md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, plan.size).gpu, content));
    md_temp_end(temp);
    md_render_mol_ctx_t kernels = md_render_mol_create(fx.dev, NULL);
    ASSERT_TRUE(kernels != NULL);
    md_gpu_addr_t dyn = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, topo.dynamic_size_velocity).gpu;

    /* the covalent bonds as the caller's pairs, their ids reversed */
    const uint32_t kv = MD_RENDER_MOL_BOND_KIND(MD_BOND_KIND_COVALENT);
    uint32_t np = 0;
    uint32_t* pairs = malloc(sys.bond.count * 8);
    for (size_t i = 0; i < sys.bond.count; ++i) {
        if (MD_RENDER_MOL_BOND_KIND(md_bond_kind(sys.bond.flags[i])) != kv) continue;
        pairs[2 * np] = sys.bond.pairs[i].idx[0]; pairs[2 * np + 1] = sys.bond.pairs[i].idx[1]; np++;
    }
    ASSERT_GT(np, 1000u);
    uint32_t* ids = malloc(np * 4);
    for (uint32_t i = 0; i < np; ++i) ids[i] = np - 1 - i;
    uint32_t offset[MD_RENDER_MOL_BOND_KIND_COUNT + 1];
    for (uint32_t k = 0; k <= MD_RENDER_MOL_BOND_KIND_COUNT; ++k) offset[k] = k <= kv ? 0 : np;

    float cell[9];
    md_render_mol_frame_cell(cell, &state.unitcell);
    float* x = malloc((size_t)n * 12);
    memcpy(x, state.xyz, (size_t)n * 12);
    md_gpu_addr_t P = f_buf(&fx, x, (size_t)n * 12);
    md_render_mol_atoms_t atoms = {0};
    atoms.position = P;
    memcpy(atoms.cell, cell, sizeof(cell));
    atoms.velocity_dt = 1.0f;

    md_render_mol_rep_t rp = {0}, rb = {0};
    rp.key = 7; rp.type = MD_RENDER_MOL_REP_PAIRS; rp.params.pairs.radius = 0.25f;
    rp.pairs = f_buf(&fx, pairs, (size_t)np * 8); rp.pair_id = f_buf(&fx, ids, (size_t)np * 4); rp.pair_offset = offset;
    rb.key = 8; rb.type = MD_RENDER_MOL_REP_BONDS; rb.params.bonds.radius = 0.25f; rb.bond_kinds = MD_RENDER_MOL_BONDS(MD_BOND_KIND_COVALENT);

    float lo[3] = { 1e30f, 1e30f, 1e30f }, hi[3] = { -1e30f, -1e30f, -1e30f };
    for (uint32_t i = 0; i < n; ++i) for (int k = 0; k < 3; ++k) { lo[k] = fminf(lo[k], x[3 * i + k]); hi[k] = fmaxf(hi[k], x[3 * i + k]); }
    const uint32_t W = 192, H = 192;
    cam_t c = { .ortho = true, .half_h = 0.6f * fmaxf(hi[0] - lo[0], hi[1] - lo[1]) + 2.0f, .aspect = 1.0f, .n = 0.0f, .far_ = 1000.0f };
    md_render_camera_t cam = make_camera(&c);
    cam.view[12] = -0.5f * (lo[0] + hi[0]); cam.view[13] = -0.5f * (lo[1] + hi[1]); cam.view[14] = -(hi[2] + 40.0f);
    md_render_view_t vb = md_render_view_create(fx.r, "bonds");
    ASSERT_TRUE(vb != NULL);
    float* dp = malloc(W * H * 4); float* db = malloc(W * H * 4);
    uint32_t* mp = malloc(W * H * 4); uint32_t* mb = malloc(W * H * 4);

    for (int step = 0; step < 2; ++step) {
        if (step == 1) {
            /* every atom a little further, every fifth small molecule across the cell too */
            for (uint32_t i = 0; i < n; ++i) { x[3 * i + 0] += 0.3f; x[3 * i + 1] -= 0.1f * (float)(i % 3); }
            const md_molecule_data_t* M = &sys.molecule;
            for (size_t m = 0; m < M->count; ++m) {
                if (M->offset[m + 1] - M->offset[m] > MD_RENDER_MOL_SMALL || m % 5) continue;
                for (uint32_t s = M->offset[m]; s < M->offset[m + 1]; ++s) for (int k = 0; k < 3; ++k) x[3 * M->atom_idx[s] + k] += cell[k];
            }
            md_gpu_upload(fx.gfx, P, x, (size_t)n * 12);
        }
        ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, MD_RENDER_MOL_UPDATE_REGROUP, content));
        md_render_batch_t bp[8], bb[8];
        const uint32_t np_b = md_render_mol_emit(&topo, dyn, &atoms, &rp, bp, 8);
        const uint32_t nb_b = md_render_mol_emit(&topo, dyn, &atoms, &rb, bb, 8);
        ASSERT_EQ(1u, np_b);
        ASSERT_EQ(1u, nb_b);
        EXPECT_EQ((uint64_t)(7 << 8 | (MD_RENDER_MOL_PART_PAIRS + kv)), bp[0].key);
        EXPECT_TRUE(bp[0].points.velocity != 0);
        md_render_frame_t fp = frame_of(W, H, cam, bp, 1), fb = frame_of(W, H, cam, bb, 1);
        memcpy(fp.content, content, sizeof(content)); memcpy(fb.content, content, sizeof(content));
        fp.time = fb.time = (double)step;
        ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &fp, NULL));
        ASSERT_TRUE(md_render_view_draw(fx.gfx, vb, &fb, NULL));
        ASSERT_TRUE(f_read(&fx, md_render_view_targets(fx.view).depth, dp, W * H * 4) && f_read(&fx, md_render_view_targets(vb).depth, db, W * H * 4));
        ASSERT_TRUE(f_read(&fx, md_render_view_targets(fx.view).motion, mp, W * H * 4) && f_read(&fx, md_render_view_targets(vb).motion, mb, W * H * 4));
        int covered = 0, depth_off = 0, motion_off = 0, moving = 0;
        for (uint32_t i = 0; i < W * H; ++i) {
            covered += db[i] > 0.0f;
            depth_off += dp[i] != db[i];
            motion_off += db[i] > 0.0f && mp[i] != mb[i];
            moving += db[i] > 0.0f && mb[i] != 0;
        }
        printf("    step %d: %d pixels of bonds, %d moving; %d of another depth, %d of another motion\n", step, covered, moving, depth_off, motion_off);
        EXPECT_GT(covered, (int)(W * H / 20));
        EXPECT_EQ(0, depth_off);
        EXPECT_EQ(0, motion_off);
        if (step == 1) EXPECT_GT(moving, covered / 2);
    }

    /* a pick on a pair: its id, the pair the topology's bond there joins */
    uint32_t px = W, py = H;
    for (uint32_t i = 0; i < W * H && px == W; ++i) if (dp[i] > 0.0f && (i % W) > W / 3 && (i / W) > H / 3) { px = i % W; py = i / W; }
    ASSERT_LT(px, W);
    md_render_pick_t pk;
    md_render_mol_hit_t hp, hb;
    ASSERT_TRUE(f_pick(&fx, px, py, &pk));
    ASSERT_TRUE(md_render_mol_decode_pick(&pk, &hp));
    EXPECT_EQ(MD_RENDER_MOL_HIT_PAIR, hp.kind);
    EXPECT_EQ(7u, (unsigned)hp.rep_key);
    ASSERT_LT(hp.bond, np);
    const uint32_t pi = np - 1 - hp.bond;                      /* the id's pair */
    EXPECT_TRUE(hp.atom == pairs[2 * pi] || hp.atom == pairs[2 * pi + 1]);
    ASSERT_TRUE(md_render_view_pick(fx.gfx, vb, px, py));
    md_gpu_stream_sync(fx.gfx);
    ASSERT_TRUE(md_render_view_pick_result(vb, &pk));
    ASSERT_TRUE(md_render_mol_decode_pick(&pk, &hb));
    EXPECT_EQ(MD_RENDER_MOL_HIT_BOND, hb.kind);
    const md_atom_pair_t bp2 = sys.bond.pairs[hb.bond];
    EXPECT_TRUE((bp2.idx[0] == pairs[2 * pi] && bp2.idx[1] == pairs[2 * pi + 1]) || (bp2.idx[1] == pairs[2 * pi] && bp2.idx[0] == pairs[2 * pi + 1]));

    /* refused: no pairs, or no offsets */
    md_render_mol_rep_t bad = rp;
    bad.pair_offset = NULL;
    EXPECT_EQ(0u, md_render_mol_emit(&topo, dyn, &atoms, &bad, NULL, 0));

    free(dp); free(db); free(mp); free(mb); free(x); free(pairs); free(ids);
    md_render_view_destroy(fx.gfx, vb);
    md_gpu_free(fx.gfx, P);
    md_gpu_free(fx.gfx, dyn);
    md_gpu_free(fx.gfx, topo.base);
    md_render_mol_destroy(kernels);
    md_vm_arena_destroy(alloc);
    f_close(&fx);
}

/* A representation's batches drawn into a view, its depth and motion read */
static bool mol_draw(fix_t* fx, md_render_view_t v, const md_render_mol_topo_t* topo, md_gpu_addr_t dyn, const md_render_mol_atoms_t* atoms,
                     const md_render_mol_rep_t* rep, md_render_camera_t cam, uint32_t W, uint32_t H, const uint64_t* content, double time,
                     float* depth, uint32_t* motion) {
    md_render_batch_t b[32];
    const uint32_t n = md_render_mol_emit(topo, dyn, atoms, rep, b, 32);
    if (n == 0 || n > 32) return false;
    md_render_frame_t f = frame_of(W, H, cam, b, n);
    memcpy(f.content, content, sizeof(f.content));
    f.time = time;
    return md_render_view_draw(fx->gfx, v, &f, NULL) && f_read(fx, md_render_view_targets(v).depth, depth, W * H * 4) &&
           f_read(fx, md_render_view_targets(v).motion, motion, W * H * 4);
}

/* Depth and motion of two draws compared: pixels of another depth, of
   another motion where the depth is the same (by more than a twentieth of a
   pixel: a tie at a joint between two sticks' round ends goes to either,
   llvmpipe from run to run, and their motions there differ by a rounding) */
static void depth_motion_off(const float* d0, const uint32_t* m0, const float* d1, const uint32_t* m1, uint32_t n, int* covered, int* depth_off,
                             int* motion_off) {
    *covered = *depth_off = *motion_off = 0;
    for (uint32_t i = 0; i < n; ++i) {
        *covered   += d0[i] > 0.0f;
        *depth_off += d1[i] != d0[i];
        if (!(d0[i] > 0.0f) || d1[i] != d0[i] || m1[i] == m0[i]) continue;
        const float dm = fmaxf(fabsf(half_bits_to_float(m1[i] & 0xFFFFu) - half_bits_to_float(m0[i] & 0xFFFFu)),
                               fabsf(half_bits_to_float(m1[i] >> 16) - half_bits_to_float(m0[i] >> 16)));
        *motion_off += dm > 0.05f;
    }
}

/* A frame's bonds (md_render_mol_bonds_t) that are the system's draw what
   the system draws, through a regroup that moves the waters (and their
   bonds). Without a chemistry: every third bond broken and formed again,
   the sticks (a pick of one formed reports its id, sys->bond.count past
   its index, never a broken bond); only broken, less. With the whole system
   as the region (md_chem_delta_t of every atom, bond and ring, as perceived)
   and the bond orders: its sticks, multiple sticks, circles and partners,
   all of them the frame's own */
UTEST(render, mol_frame_bonds_draw_the_system) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(load_system(&sys, &state, MD_UNITTEST_DATA_DIR "/tryptophan-md.gro"));
    const uint32_t n = (uint32_t)sys.atom.count, count = (uint32_t)sys.bond.count;
    uint64_t content[MD_RENDER_CONTENT_COUNT] = {0};
    md_temp_scope_t temp = md_temp_begin();
    md_render_mol_plan_t plan;
    ASSERT_TRUE(md_render_mol_plan(&plan, &sys, NULL, temp));
    md_render_mol_topo_t topo;
    ASSERT_TRUE(md_render_mol_write(fx.gfx, &topo, &plan, md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, plan.size).gpu, content));
    ASSERT_GT(topo.small_bond_count, 100u);
    ASSERT_GT(topo.orders.order_count[0], 0u);
    ASSERT_GT(topo.orders.aromatic_count[3] + topo.orders.aromatic_count[2], 0u);
    md_render_mol_ctx_t kernels = md_render_mol_create(fx.dev, NULL);
    ASSERT_TRUE(kernels != NULL);
    md_gpu_addr_t dyn = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, topo.dynamic_size_velocity).gpu;

    /* every third covalent bond broken and formed again (its flags kept) */
    md_bitmask_t broken = md_bitmask_create(0, count, alloc);
    md_atom_pair_t*  fpair = md_alloc(alloc, ((size_t)count + 1) * sizeof(md_atom_pair_t));
    md_bond_flags_t* fflag = md_alloc(alloc, ((size_t)count + 1) * sizeof(md_bond_flags_t));
    uint32_t nf = 0;
    for (uint32_t i = 0; i < count; ++i) {
        if (i % 3 || md_bond_kind(sys.bond.flags[i]) != MD_BOND_KIND_COVALENT) continue;
        md_bitmask_set_bit(&broken, i);
        fpair[nf] = sys.bond.pairs[i];
        fflag[nf++] = sys.bond.flags[i];
    }
    ASSERT_GT(nf, 1000u);
    md_bond_delta_t again  = { .broken = md_bitmask_view(&broken), .num_formed = nf, .formed_pair = fpair, .formed_flags = fflag };
    md_bond_delta_t gone   = { .broken = md_bitmask_view(&broken) };
    /* the whole system as the region, as perceived */
    md_chem_delta_t whole = {0};
    whole.num_atoms  = n;
    whole.atom_idx   = md_alloc(alloc, (size_t)n * sizeof(int32_t));
    whole.atom_flags = md_alloc(alloc, (size_t)n * sizeof(md_atom_flags_t));
    for (uint32_t i = 0; i < n; ++i) { whole.atom_idx[i] = (int32_t)i; whole.atom_flags[i] = sys.atom.flags[i]; }
    whole.bond_idx   = md_alloc(alloc, ((size_t)count + 1) * 4u);
    whole.bond_flags = md_alloc(alloc, ((size_t)count + 1) * sizeof(md_bond_flags_t));
    for (uint32_t i = 0; i < count; ++i) {
        if (md_bond_kind(sys.bond.flags[i]) != MD_BOND_KIND_COVALENT) continue;
        whole.bond_idx[whole.num_bonds] = i;
        whole.bond_flags[whole.num_bonds++] = sys.bond.flags[i];
    }
    whole.ring = sys.ring;
    const md_bond_delta_t none = {0};

    md_render_mol_bonds_t fb[3];
    const md_bond_delta_t* deltas[3] = { &again, &gone, &none };
    const md_chem_delta_t* chems[3]  = { NULL, NULL, &whole };
    for (int k = 0; k < 3; ++k) {
        md_render_mol_bonds_plan_t bp;
        ASSERT_TRUE(md_render_mol_bonds_plan(&bp, &sys, NULL, NULL, deltas[k], chems[k], temp));
        ASSERT_TRUE(md_render_mol_bonds_write(fx.gfx, &fb[k], &bp, md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, bp.size).gpu, content));
    }
    EXPECT_EQ(nf, fb[0].pair_offset[MD_RENDER_MOL_BOND_KIND_COUNT]);
    EXPECT_EQ(0u, fb[1].pair_offset[MD_RENDER_MOL_BOND_KIND_COUNT]);
    EXPECT_EQ(topo.orders.order_count[0], fb[2].orders.order_count[0]);
    EXPECT_EQ(topo.orders.order_count[1], fb[2].orders.order_count[1]);
    EXPECT_EQ(topo.orders.partner_count, fb[2].orders.partner_count);
    for (uint32_t k = 0; k < MD_RENDER_MOL_RING_SIZES; ++k) EXPECT_EQ(topo.orders.aromatic_count[k], fb[2].orders.aromatic_count[k]);
    md_temp_end(temp);

    float cell[9];
    md_render_mol_frame_cell(cell, &state.unitcell);
    float* x = malloc((size_t)n * 12);
    memcpy(x, state.xyz, (size_t)n * 12);
    md_gpu_addr_t P = f_buf(&fx, x, (size_t)n * 12);
    md_render_mol_atoms_t atoms = {0};
    atoms.position = P;
    memcpy(atoms.cell, cell, sizeof(cell));
    atoms.velocity_dt = 1.0f;

    md_render_mol_rep_t rep = {0};
    rep.key = 5; rep.type = MD_RENDER_MOL_REP_LICORICE; rep.params.licorice.radius = 0.25f;
    float lo[3] = { 1e30f, 1e30f, 1e30f }, hi[3] = { -1e30f, -1e30f, -1e30f };
    for (uint32_t i = 0; i < n; ++i) for (int k = 0; k < 3; ++k) { lo[k] = fminf(lo[k], x[3 * i + k]); hi[k] = fmaxf(hi[k], x[3 * i + k]); }
    const uint32_t W = 192, H = 192;
    cam_t c = { .ortho = true, .half_h = 0.6f * fmaxf(hi[0] - lo[0], hi[1] - lo[1]) + 2.0f, .aspect = 1.0f, .n = 0.0f, .far_ = 1000.0f };
    md_render_camera_t cam = make_camera(&c);
    cam.view[12] = -0.5f * (lo[0] + hi[0]); cam.view[13] = -0.5f * (lo[1] + hi[1]); cam.view[14] = -(hi[2] + 40.0f);
    /* a view per draw, each drawn once a step: the motion over the same time */
    md_render_view_t vs[2], vf[2][3];
    for (int o = 0; o < 2; ++o) {
        vs[o] = md_render_view_create(fx.r, "system");
        ASSERT_TRUE(vs[o] != NULL);
        for (int k = 0; k < 3; ++k) { vf[o][k] = md_render_view_create(fx.r, "frame"); ASSERT_TRUE(vf[o][k] != NULL); }
    }
    float* d0 = malloc(W * H * 4); float* d1 = malloc(W * H * 4);
    uint32_t* m0 = malloc(W * H * 4); uint32_t* m1 = malloc(W * H * 4);

    /* a frame, then every atom a little further and every fifth small
       molecule across the cell: regrouped */
    for (int step = 0; step < 2; ++step) {
        if (step == 1) {
            for (uint32_t i = 0; i < n; ++i) { x[3 * i + 0] += 0.3f; x[3 * i + 1] -= 0.1f * (float)(i % 3); }
            const md_molecule_data_t* M = &sys.molecule;
            for (size_t m = 0; m < M->count; ++m) {
                if (M->offset[m + 1] - M->offset[m] > MD_RENDER_MOL_SMALL || m % 5) continue;
                for (uint32_t q = M->offset[m]; q < M->offset[m + 1]; ++q) for (int k = 0; k < 3; ++k) x[3 * M->atom_idx[q] + k] += cell[k];
            }
            md_gpu_upload(fx.gfx, P, x, (size_t)n * 12);
        }
        /* one update a step, with the bonds whose orders want sides (theirs
           kept in their block, so the other draws need none) */
        atoms.bonds = &fb[2];
        ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, MD_RENDER_MOL_UPDATE_REGROUP, content));
        for (int orders = 0; orders < 2; ++orders) {
            rep.bond_orders = orders != 0;
            atoms.bonds = NULL;
            ASSERT_TRUE(mol_draw(&fx, vs[orders], &topo, dyn, &atoms, &rep, cam, W, H, content, step, d0, m0));
            int moving = 0;
            for (uint32_t i = 0; i < W * H; ++i) moving += d0[i] > 0.0f && m0[i] != 0;
            if (step == 1) EXPECT_GT(moving, (int)(W * H / 40));
            for (int k = 0; k < 3; ++k) {
                if (orders != (k == 2) && k != 1) continue;
                atoms.bonds = &fb[k];
                ASSERT_TRUE(mol_draw(&fx, vf[orders][k], &topo, dyn, &atoms, &rep, cam, W, H, content, step, d1, m1));
                int covered, depth_off, motion_off;
                depth_motion_off(d0, m0, d1, m1, W * H, &covered, &depth_off, &motion_off);
                static const char* what[3] = { "every third bond broken and formed again", "every third bond broken", "the region the whole system" };
                printf("    step %d, %s, %s: %d of %d pixels of another depth (%d moving), %d of another motion\n", step, orders ? "bond orders" : "sticks",
                       what[k], depth_off, covered, moving, motion_off);
                EXPECT_GT(covered, (int)(W * H / 20));
                if (k == 1) {
                    EXPECT_GT(depth_off, covered / 50);
                    continue;
                }
                EXPECT_EQ(0, depth_off);
                EXPECT_EQ(0, motion_off);
                if (k != 0 || step != 1) continue;
                /* picks: a formed bond reports its id, never a broken one */
                int formed_hits = 0, broken_hits = 0, bond_hits = 0;
                for (uint32_t i = 0; i < W * H && formed_hits < 4; i += 37) {
                    if (!(d1[i] > 0.0f)) continue;
                    md_render_pick_t pk;
                    md_render_mol_hit_t h;
                    if (!md_render_view_pick(fx.gfx, vf[orders][k], i % W, i / W)) continue;
                    md_gpu_stream_sync(fx.gfx);
                    if (!md_render_view_pick_result(vf[orders][k], &pk) || !md_render_mol_decode_pick(&pk, &h)) continue;
                    if (h.kind == MD_RENDER_MOL_HIT_PAIR) {
                        formed_hits++;
                        ASSERT_GE(h.bond, count);
                        ASSERT_LT(h.bond - count, nf);
                        const md_atom_pair_t bp = fpair[h.bond - count];
                        EXPECT_TRUE(h.atom == (uint32_t)bp.idx[0] || h.atom == (uint32_t)bp.idx[1]);
                    } else if (h.kind == MD_RENDER_MOL_HIT_BOND) {
                        bond_hits++;
                        broken_hits += h.bond < count && md_bitmask_test_bit(&broken, h.bond);
                    }
                }
                printf("    picks: %d of formed bonds, %d of bonds, %d of broken ones\n", formed_hits, bond_hits, broken_hits);
                EXPECT_GT(formed_hits, 0);
                EXPECT_GT(bond_hits, 0);
                EXPECT_EQ(0, broken_hits);
            }
        }
    }
    /* refused: pairs on anything but PAIRS */
    md_render_mol_rep_t bad = rep;
    uint32_t offs[MD_RENDER_MOL_BOND_KIND_COUNT + 1] = {0};
    bad.pairs = fb[0].pairs; bad.pair_offset = offs;
    EXPECT_EQ(0u, md_render_mol_emit(&topo, dyn, &atoms, &bad, NULL, 0));

    free(d0); free(d1); free(m0); free(m1); free(x);
    for (int o = 0; o < 2; ++o) {
        md_render_view_destroy(fx.gfx, vs[o]);
        for (int k = 0; k < 3; ++k) md_render_view_destroy(fx.gfx, vf[o][k]);
    }
    md_gpu_free(fx.gfx, P);
    for (int k = 0; k < 3; ++k) md_gpu_free(fx.gfx, fb[k].bond_visible);
    md_gpu_free(fx.gfx, dyn);
    md_gpu_free(fx.gfx, topo.base);
    md_render_mol_destroy(kernels);
    md_vm_arena_destroy(alloc);
    f_close(&fx);
}

/* Benzene, a small molecule (12 atoms): its circle drawn (md_render_mol's
   small molecules regroup, their orders through each atom's point), picked
   where the molecule is, before and after the regroup moves them */
UTEST(render, mol_small_molecule_orders) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    /* 27 benzenes in the xy plane, 9 apart in a box of 27 */
    char* pdb = md_alloc(alloc, 1 << 16);
    int len = sprintf(pdb, "CRYST1   27.000   27.000   27.000  90.00  90.00  90.00 P 1           1\n");
    int serial = 1;
    for (int m = 0; m < 27; ++m) {
        const float cx = 4.5f + 9.0f * (float)(m % 3), cy = 4.5f + 9.0f * (float)((m / 3) % 3), cz = 4.5f + 9.0f * (float)(m / 9);
        for (int e = 0; e < 2; ++e) {
            for (int k = 0; k < 6; ++k) {
                const float a = 1.0471976f * (float)k, r = e ? 2.479f : 1.39f;
                len += sprintf(pdb + len, "HETATM%5d %-4s BNZ A%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n", serial++, e ? " H" : " C", m + 1,
                               cx + r * cosf(a), cy + r * sinf(a), cz, e ? "H" : "C");
            }
        }
    }
    len += sprintf(pdb + len, "END\n");
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, (md_str_t){ pdb, (size_t)len }, NULL));
    ASSERT_TRUE(md_system_infer(&sys, &state, NULL));
    const uint32_t n = (uint32_t)sys.atom.count;
    ASSERT_EQ(27u * 12u, n);
    uint64_t content[MD_RENDER_CONTENT_COUNT] = {0};
    md_temp_scope_t temp = md_temp_begin();
    md_render_mol_plan_t plan;
    ASSERT_TRUE(md_render_mol_plan(&plan, &sys, NULL, temp));
    md_render_mol_topo_t topo;
    ASSERT_TRUE(md_render_mol_write(fx.gfx, &topo, &plan, md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, plan.size).gpu, content));
    md_temp_end(temp);
    EXPECT_EQ(27u, topo.small_molecule_count);
    EXPECT_EQ(27u, topo.orders.aromatic_count[3]);
    md_render_mol_ctx_t kernels = md_render_mol_create(fx.dev, NULL);
    ASSERT_TRUE(kernels != NULL);
    md_gpu_addr_t dyn = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, topo.dynamic_size).gpu;
    float* x = malloc((size_t)n * 12);
    memcpy(x, state.xyz, (size_t)n * 12);
    md_gpu_addr_t P = f_buf(&fx, x, (size_t)n * 12);
    md_render_mol_atoms_t atoms = {0};
    atoms.position = P;
    md_render_mol_frame_cell(atoms.cell, &state.unitcell);
    md_render_mol_rep_t rep = {0};
    rep.key = 3; rep.type = MD_RENDER_MOL_REP_LICORICE; rep.params.licorice.radius = 0.2f; rep.bond_orders = true;
    const uint32_t W = 256, H = 256;
    cam_t c = { .ortho = true, .half_h = 16.0f, .aspect = 1.0f, .n = 0.0f, .far_ = 1000.0f };
    md_render_camera_t cam = make_camera(&c);
    cam.view[12] = -13.5f; cam.view[13] = -13.5f; cam.view[14] = -100.0f;
    float* depth = malloc(W * H * 4);
    uint32_t* motion = malloc(W * H * 4);
    /* each benzene of the front layer (z 22.5) picked on its circle: the
       apothem 1.204 less the inset 0.38, towards an edge's middle */
    for (int step = 0; step < 2; ++step) {
        if (step == 1) {        /* every molecule a third of the box along x, wrapped: regrouped */
            for (uint32_t i = 0; i < n; ++i) { x[3 * i] += 9.0f; if (x[3 * i] > 27.0f) x[3 * i] -= 27.0f; }
            md_gpu_upload(fx.gfx, P, x, (size_t)n * 12);
        }
        ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, MD_RENDER_MOL_UPDATE_REGROUP, content));
        ASSERT_TRUE(mol_draw(&fx, fx.view, &topo, dyn, &atoms, &rep, cam, W, H, content, step, depth, motion));
        int rings = 0, right = 0;
        for (int m = 18; m < 27; ++m) {
            float cx = 0, cy = 0;
            for (int k = 0; k < 6; ++k) { cx += x[3 * (12 * m + k)] / 6.0f; cy += x[3 * (12 * m + k) + 1] / 6.0f; }
            const float px = cx + 0.824f * cosf(0.5235988f), py = cy + 0.824f * sinf(0.5235988f);
            md_render_pick_t pk;
            md_render_mol_hit_t h;
            const uint32_t sx = (uint32_t)(((px + cam.view[12]) / c.half_h * 0.5f + 0.5f) * W);
            const uint32_t sy = (uint32_t)((0.5f - (py + cam.view[13]) / c.half_h * 0.5f) * H);
            if (!f_pick(&fx, sx, sy, &pk) || !md_render_mol_decode_pick(&pk, &h)) continue;
            rings += h.kind == MD_RENDER_MOL_HIT_RING;
            right += h.kind == MD_RENDER_MOL_HIT_RING && h.atom / 12u == (uint32_t)m && md_ring_atoms(&sys.ring, h.ring)[0] / 12 == m;
        }
        printf("    step %d: %d of the 9 front benzenes' circles picked, %d as their own\n", step, rings, right);
        EXPECT_EQ(9, rings);
        EXPECT_EQ(9, right);
    }
    free(x); free(depth); free(motion);
    md_gpu_free(fx.gfx, P);
    md_gpu_free(fx.gfx, dyn);
    md_gpu_free(fx.gfx, topo.base);
    md_render_mol_destroy(kernels);
    md_vm_arena_destroy(alloc);
    f_close(&fx);
}

/* A Diels-Alder reaction (test_chem's: butadiene and ethylene, as a frame
   cyclohexene): the frame's bonds as md_chem perceives them. The formed
   bonds and the region's are the frame's own: one double bond, the middle
   one of the butadiene, which the system has single; the system's three
   doubles hidden. Drawn, its sticks pick as that bond */
static const char render_diels_alder_before[] =
    "16\n" "diels_alder before\n"
    "C 0.8500 -0.0113 1.0574\n" "C 1.9229 0.4135 0.3803\n" "C 3.0805 -0.4135 0.1454\n" "C 4.1534 0.0113 -0.5317\n"
    "C 11.9048 0.0070 0.0359\n" "C 10.5710 -0.0070 -0.0359\n" "H 0.7872 -1.0117 1.4734\n" "H 0.0000 0.6463 1.2104\n"
    "H 1.9319 1.4272 -0.0150\n" "H 3.0715 -1.4272 0.5408\n" "H 5.0034 -0.6463 -0.6847\n" "H 4.2162 1.0117 -0.9476\n"
    "H 12.4725 -0.9167 0.0885\n" "H 12.4553 0.9426 0.0434\n" "H 10.0206 -0.9426 -0.0434\n" "H 10.0034 0.9167 -0.0885\n";
static const float render_diels_alder_after[16][3] = {
    { 3.4709f, -0.6842f, 0.2563f }, { 3.5206f, 0.6899f, -0.3459f }, { 2.4234f, 1.4151f, -0.6089f }, { 1.0302f, 0.9144f, -0.3607f },
    { 0.9723f, -0.5923f, -0.0956f }, { 2.0979f, -1.0368f, 0.8345f }, { 4.2226f, -0.7617f, 1.0500f }, { 3.7444f, -1.4108f, -0.5180f },
    { 4.5024f, 1.1004f, -0.5682f }, { 2.5261f, 2.4139f, -1.0257f }, { 0.3992f, 1.1488f, -1.2256f }, { 0.6155f, 1.4576f, 0.4968f },
    { 1.0530f, -1.1351f, -1.0462f }, { 0.0000f, -0.8534f, 0.3374f }, { 2.0379f, -2.1168f, 1.0098f }, { 1.9748f, -0.5491f, 1.8101f } };

static int32_t find_bond_of(const md_system_t* sys, int a, int b) {
    for (size_t i = 0; i < sys->bond.count; ++i) {
        const md_atom_pair_t p = sys->bond.pairs[i];
        if ((p.idx[0] == a && p.idx[1] == b) || (p.idx[0] == b && p.idx[1] == a)) return (int32_t)i;
    }
    return -1;
}

UTEST(render, mol_frame_bonds_diels_alder) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_xyz_system_init_from_str(&sys, &state, md_str_from_cstr(render_diels_alder_before), NULL));
    ASSERT_TRUE(md_system_infer(&sys, &state, NULL));
    md_system_state_t frame = { .alloc = alloc };
    md_system_state_init(&frame, 16);
    for (int i = 0; i < 16; ++i) frame.xyz[i] = (md_float3_t){ render_diels_alder_after[i][0], render_diels_alder_after[i][1], render_diels_alder_after[i][2] };
    md_bond_delta_t delta = {0};
    const md_bond_delta_desc_t desc = { .tolerance = MD_BOND_DELTA_TOLERANCE, .breakable = MD_BOND_DELTA_BREAKABLE };
    ASSERT_TRUE(md_bond_delta_compute(&delta, &sys, &frame, &desc, alloc));
    ASSERT_EQ(2, (int)delta.num_formed);
    md_chem_delta_t chem = {0};
    ASSERT_TRUE(md_chem_perceive_delta(&chem, &sys, &frame, &delta, MD_CHEM_FLAG_PROTONATE_PH7, alloc));
    const int32_t c12 = find_bond_of(&sys, 1, 2), c01 = find_bond_of(&sys, 0, 1), c23 = find_bond_of(&sys, 2, 3), c45 = find_bond_of(&sys, 4, 5);
    ASSERT_TRUE(c12 >= 0 && c01 >= 0 && c23 >= 0 && c45 >= 0);
    const uint32_t count = (uint32_t)sys.bond.count;

    md_temp_scope_t temp = md_temp_begin();
    md_render_mol_plan_t plan;
    ASSERT_TRUE(md_render_mol_plan(&plan, &sys, NULL, temp));
    EXPECT_EQ(3u, plan.layout.orders.order_count[0]);       /* the system's: butadiene's two, ethylene's */
    md_render_mol_bonds_plan_t bp;
    ASSERT_TRUE(md_render_mol_bonds_plan(&bp, &sys, NULL, &frame, &delta, &chem, temp));
    const md_render_mol_bonds_t* L = &bp.layout;
    const uint8_t* img = bp.data;
    ASSERT_EQ(1u, L->orders.order_count[0]);
    EXPECT_EQ(0u, L->orders.order_count[1]);
    EXPECT_EQ((uint32_t)c12, ((const uint32_t*)(img + L->orders.order_id))[0]);
    const uint32_t* oa = (const uint32_t*)(img + L->orders.order_atoms);
    EXPECT_TRUE((oa[0] == 1 && oa[1] == 2) || (oa[0] == 2 && oa[1] == 1));
    const uint32_t* bv = (const uint32_t*)(img + L->bond_visible);
    #define BIT(w, i) (((w)[(i) >> 5] >> ((i) & 31)) & 1u)
    EXPECT_EQ(0u, BIT(bv, (uint32_t)c12));
    EXPECT_EQ(0u, BIT(bv, (uint32_t)c01));
    EXPECT_EQ(0u, BIT(bv, (uint32_t)c23));
    EXPECT_EQ(0u, BIT(bv, (uint32_t)c45));
    #undef BIT
    /* the frame's own sticks: every bond of the region but the double, the formed ones among them by their ids */
    const uint32_t np = L->pair_offset[MD_RENDER_MOL_BOND_KIND_COUNT];
    const uint32_t* pid = (const uint32_t*)(img + L->pair_id);
    uint32_t formed_ids = 0;
    for (uint32_t k = 0; k < np; ++k) { formed_ids += pid[k] >= count; EXPECT_NE((uint32_t)c12, pid[k]); }
    printf("    the frame's own: %u sticks (%u formed), %u double\n", np, formed_ids, L->orders.order_count[0]);
    EXPECT_EQ(2u, formed_ids);
    EXPECT_EQ((uint32_t)chem.num_bonds - 1u, np);

    fix_t fx;
    if (!f_open(&fx)) { md_temp_end(temp); md_vm_arena_destroy(alloc); UTEST_SKIP(f_skip_reason()); }
    uint64_t content[MD_RENDER_CONTENT_COUNT] = {0};
    md_render_mol_topo_t topo;
    ASSERT_TRUE(md_render_mol_write(fx.gfx, &topo, &plan, md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, plan.size).gpu, content));
    md_render_mol_bonds_t fb;
    ASSERT_TRUE(md_render_mol_bonds_write(fx.gfx, &fb, &bp, md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, bp.size).gpu, content));
    md_temp_end(temp);
    md_render_mol_ctx_t kernels = md_render_mol_create(fx.dev, NULL);
    ASSERT_TRUE(kernels != NULL);
    md_gpu_addr_t dyn = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, topo.dynamic_size).gpu;
    md_gpu_addr_t P = f_buf(&fx, frame.xyz, 16 * 12);
    md_render_mol_atoms_t atoms = {0};
    atoms.position = P;
    atoms.bonds = &fb;
    ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, MD_RENDER_MOL_UPDATE_REGROUP, content));
    md_render_mol_rep_t rep = {0};
    rep.key = 9; rep.type = MD_RENDER_MOL_REP_LICORICE; rep.params.licorice.radius = 0.2f; rep.bond_orders = true;
    /* looking down the plane of the ring's double bond face on: along its
       neighbours' normal, so its two sticks lie side by side */
    const md_float3_t a = frame.xyz[1], b = frame.xyz[2], r3 = frame.xyz[3];
    const double u[3] = { b.x - a.x, b.y - a.y, b.z - a.z }, v[3] = { r3.x - a.x, r3.y - a.y, r3.z - a.z };
    double nz[3];
    cross3(u, v, nz);
    double ln = sqrt(dot3(nz, nz)), lu = sqrt(dot3(u, u));
    double ex[3] = { u[0] / lu, u[1] / lu, u[2] / lu }, ez[3] = { nz[0] / ln, nz[1] / ln, nz[2] / ln }, ey[3];
    cross3(ez, ex, ey);
    const double mid[3] = { 0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z) };
    const uint32_t W = 128, H = 128;
    cam_t c = { .ortho = true, .half_h = 2.0f, .aspect = 1.0f, .n = 0.0f, .far_ = 1000.0f };
    md_render_camera_t cam = make_camera(&c);
    for (int i = 0; i < 3; ++i) { cam.view[4 * i + 0] = (float)ex[i]; cam.view[4 * i + 1] = (float)ey[i]; cam.view[4 * i + 2] = (float)ez[i]; }
    cam.view[12] = (float)-dot3(ex, mid); cam.view[13] = (float)-dot3(ey, mid); cam.view[14] = (float)(-dot3(ez, mid) - 50.0);
    float* depth = malloc(W * H * 4);
    uint32_t* motion = malloc(W * H * 4);
    ASSERT_TRUE(mol_draw(&fx, fx.view, &topo, dyn, &atoms, &rep, cam, W, H, content, 0.0, depth, motion));
    /* along a vertical line through the bond's middle: the double's two sticks, as the bond */
    int on = 0, hidden = 0;
    for (uint32_t py = 0; py < H; ++py) {
        md_render_pick_t pk;
        md_render_mol_hit_t h;
        if (!f_pick(&fx, W / 2, py, &pk) || !md_render_mol_decode_pick(&pk, &h)) continue;
        on += h.kind == MD_RENDER_MOL_HIT_PAIR && h.bond == (uint32_t)c12;
        hidden += h.kind == MD_RENDER_MOL_HIT_BOND && (h.bond == (uint32_t)c12 || h.bond == (uint32_t)c01 || h.bond == (uint32_t)c23);
    }
    printf("    picks across the double bond: %d of it, %d of the system's bonds of the region\n", on, hidden);
    EXPECT_GT(on, 4);
    EXPECT_EQ(0, hidden);
    free(depth); free(motion);
    md_gpu_free(fx.gfx, P);
    md_gpu_free(fx.gfx, dyn);
    md_gpu_free(fx.gfx, fb.bond_visible);
    md_gpu_free(fx.gfx, topo.base);
    md_render_mol_destroy(kernels);
    md_vm_arena_destroy(alloc);
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
    md_render_mol_ctx_t kernels = md_render_mol_create(fx.dev, NULL);
    ASSERT_TRUE(kernels != NULL);
    md_gpu_addr_t dyn = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, topo.dynamic_size).gpu;
    md_gpu_addr_t P = f_buf(&fx, state.xyz, (size_t)n * 12);
    md_render_mol_atoms_t atoms = {0};
    atoms.position = P;
    ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, MD_RENDER_MOL_UPDATE_REGROUP, NULL));

    md_render_mol_rep_t rep = {0};
    rep.key  = 7;
    rep.type = MD_RENDER_MOL_REP_LICORICE;
    rep.params.licorice.radius = 0.3f;
    rep.bond_kinds = MD_RENDER_MOL_BONDS_ALL;
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
    md_render_debug_view(ref, MD_RENDER_DEBUG_NO_OCCLUSION | MD_RENDER_DEBUG_NO_COMPUTE_RASTER, 0.0f);
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
    md_render_mol_destroy(kernels);
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
    md_render_mol_ctx_t kernels = md_render_mol_create(fx.dev, NULL);
    ASSERT_TRUE(kernels != NULL);
    md_gpu_addr_t dyn = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, topo.dynamic_size).gpu;
    ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, 0, content));
    ASSERT_TRUE(md_system_state_backbone_compute(&state, &sys));
    const md_secondary_structure_t* ss = md_system_state_secondary_structure(&state, &sys);
    ASSERT_TRUE(ss != NULL);
    md_gpu_addr_t SS = f_buf(&fx, ss, PB->segment.count * sizeof(int32_t));
    md_gpu_addr_t prof = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, md_render_mol_cartoon_size(&topo)).gpu;
    const float scale[3] = { 1.0f, 1.0f, 1.0f };
    ASSERT_TRUE(md_render_mol_cartoon_profile(kernels, fx.gfx, &topo, prof, SS, scale, content));

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
            ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, MD_RENDER_MOL_UPDATE_RESET_BACKBONE, content));
            free(last);
            last = (float*)read_block(&fx, axes, (uint64_t)ncp * 16);
            for (int step = 0; step < 40; ++step) {
                for (uint32_t i = 0; i < na; ++i) {
                    xyz[i].x = state.xyz[i].x + 0.15f * lcg_gauss(&rng);
                    xyz[i].y = state.xyz[i].y + 0.15f * lcg_gauss(&rng);
                    xyz[i].z = state.xyz[i].z + 0.15f * lcg_gauss(&rng);
                }
                md_gpu_upload(fx.gfx, P, xyz, (size_t)na * 12);
                ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, mode ? MD_RENDER_MOL_UPDATE_RESET_BACKBONE : 0, content));
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
    md_render_mol_destroy(kernels);
    md_system_free(&sys);
    f_close(&fx);
}

/* The cartoon of a nucleic acid (a nucleosome: 294 nucleotides of DNA in two
   strands, and its histones): the band through the C4', the rings of each
   nucleotide as polygons from their colour's corner, the sticks; the band's
   axes towards the bases; picks of a ring and of the band */
UTEST(render, mol_cartoon_nucleic) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(load_system(&sys, &state, MD_UNITTEST_DATA_DIR "/nucl-dna.gro"));
    const md_nucleic_backbone_data_t* NB = &sys.nucleic_backbone;
    const uint32_t nseg = (uint32_t)NB->segment.count;
    ASSERT_GT(nseg, 200u);
    uint32_t purines = 0;
    for (uint32_t g = 0; g < nseg; ++g) {
        md_str_t rn = md_system_residue_name(&sys, (size_t)NB->segment.res_idx[g]);
        purines += (md_str_contains_char(rn, 'A') || md_str_contains_char(rn, 'G')) ? 1u : 0u;
    }
    uint64_t content[MD_RENDER_CONTENT_COUNT] = {0};

    md_temp_scope_t temp = md_temp_begin();
    md_render_mol_plan_t plan;
    ASSERT_TRUE(md_render_mol_plan(&plan, &sys, NULL, temp));
    const md_render_mol_topo_t* L = &plan.layout;
    const uint8_t*  img   = (const uint8_t*)plan.data;
    const uint32_t* order = (const uint32_t*)(img + L->atom_order);
    printf("    %u nucleotides (%u purines) in %u strands: band %u, rings of 5 %u, of 6 %u, sticks %u\n", nseg, purines,
           (uint32_t)NB->range.count, L->nucleic_count, L->ring_count[2], L->ring_count[3], L->stick_count);
    /* Every nucleotide: a control point, a sugar, a ring of six, a purine's
       ring of five, a stick; the pieces the strands */
    EXPECT_EQ(nseg, L->nucleic_count);
    EXPECT_EQ(nseg + purines, L->ring_count[2]);
    EXPECT_EQ(nseg, L->ring_count[3]);
    EXPECT_EQ(nseg, L->stick_count);
    EXPECT_EQ(0u, L->ring_count[0] + L->ring_count[1] + L->ring_count[4]);
    {
        const uint32_t* idx  = (const uint32_t*)(img + L->nucleic_index);
        const uint32_t* seg  = (const uint32_t*)(img + L->nucleic_segment);
        const uint32_t* face = (const uint32_t*)(img + L->nucleic_face);
        int bad = 0, pieces = 0;
        for (uint32_t j = 0; j < L->nucleic_count; ++j) {
            const md_nucleic_acid_atoms_t at = NB->segment.atoms[seg[j]];
            if (order[idx[j] & 0x3FFFFFFFu] != (uint32_t)at.c4) bad++;
            uint32_t f = order[face[j]];        /* the glycosidic N: bonded to the C1', a nitrogen */
            if (md_bond_find(&sys.bond, (md_atom_idx_t)f, at.c1) < 0 || md_atom_atomic_number(&sys.atom, f) != 7) bad++;
            pieces += (idx[j] & MD_RENDER_SPLINE_FIRST) ? 1 : 0;
        }
        EXPECT_EQ(0, bad);
        EXPECT_EQ((int)NB->range.count, pieces);
        /* the rings: around, from the C1' of a sugar or a carbon of a base */
        for (uint32_t k = 2; k <= 3; ++k) {
            const uint32_t  nc  = 3 + k;
            const uint32_t* pts = (const uint32_t*)(img + L->ring_points[k]);
            const uint32_t* rs  = (const uint32_t*)(img + L->ring_segment[k]);
            int rbad = 0, sugars = 0;
            for (uint32_t r = 0; r < L->ring_count[k]; ++r) {
                const md_nucleic_acid_atoms_t at = NB->segment.atoms[rs[r]];
                for (uint32_t l = 0; l < nc; ++l) {
                    if (md_bond_find(&sys.bond, (md_atom_idx_t)order[pts[nc * r + l]], (md_atom_idx_t)order[pts[nc * r + (l + 1) % nc]]) < 0) rbad++;
                }
                uint32_t a0 = order[pts[nc * r]];
                bool sugar = false;
                for (uint32_t l = 0; l < nc; ++l) sugar |= order[pts[nc * r + l]] == (uint32_t)at.o4;
                sugars += sugar;
                if (sugar ? a0 != (uint32_t)at.c1 : md_atom_atomic_number(&sys.atom, a0) != 6) rbad++;
            }
            EXPECT_EQ(0, rbad);
            EXPECT_EQ(k == 2 ? (int)nseg : 0, sugars);
        }
        const uint32_t* sp = (const uint32_t*)(img + L->stick_points);
        const uint32_t* ss = (const uint32_t*)(img + L->stick_segment);
        int sbad = 0;
        for (uint32_t i = 0; i < L->stick_count; ++i) {
            if (order[sp[2 * i]] != (uint32_t)NB->segment.atoms[ss[i]].c1) sbad++;
            if (md_atom_atomic_number(&sys.atom, order[sp[2 * i + 1]]) != 7) sbad++;
        }
        EXPECT_EQ(0, sbad);
    }
    md_gpu_mem_t mem = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, plan.size);
    md_render_mol_topo_t topo;
    ASSERT_TRUE(md_render_mol_write(fx.gfx, &topo, &plan, mem.gpu, content));
    md_temp_end(temp);

    const uint32_t na = (uint32_t)sys.atom.count;
    md_gpu_addr_t P = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, (uint64_t)na * 12).gpu;
    ASSERT_TRUE(md_render_mol_upload_state(fx.gfx, P, &state, content));
    md_render_mol_atoms_t atoms = {0};
    atoms.position = P;
    md_render_mol_ctx_t kernels = md_render_mol_create(fx.dev, NULL);
    ASSERT_TRUE(kernels != NULL);
    md_gpu_addr_t dyn = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, topo.dynamic_size).gpu;
    ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, 0, content));

    /* The batches: the protein's band, the nucleic band, rings of 5 and 6, sticks */
    md_render_mol_rep_t rep = {0};
    rep.key  = 9;
    rep.type = MD_RENDER_MOL_REP_CARTOON;
    md_render_batch_t bt[8];
    ASSERT_EQ(5u, md_render_mol_emit(&topo, dyn, &atoms, &rep, bt, 8));
    EXPECT_EQ(MD_RENDER_SPLINES, bt[1].kind);
    EXPECT_TRUE(bt[1].style.splines.through);
    EXPECT_EQ(MD_RENDER_POLYGONS, bt[2].kind);
    EXPECT_EQ(5u, bt[2].style.polygons.corners);
    EXPECT_EQ(6u, bt[3].style.polygons.corners);
    EXPECT_EQ(MD_RENDER_CAPSULES, bt[4].kind);

    /* The axes: unit, across the chord, across the way to the base (so the
       faces face it) */
    {
        const uint32_t n = topo.nucleic_count;
        float*    ax   = (float*)read_block(&fx, bt[1].style.splines.axis, (uint64_t)n * 16);
        uint32_t* idx  = (uint32_t*)read_block(&fx, topo.nucleic_index, (uint64_t)n * 4);
        uint32_t* face = (uint32_t*)read_block(&fx, topo.nucleic_face, (uint64_t)n * 4);
        uint32_t* ord  = (uint32_t*)read_block(&fx, topo.atom_order, (uint64_t)topo.point_count * 4);
        int bad = 0;
        float worst = 0;
        for (uint32_t j = 0; j < n; ++j) {
            const float* a = ax + 4 * j;
            uint32_t jp = (idx[j] & MD_RENDER_SPLINE_FIRST) ? j : j - 1, jn = (idx[j] & MD_RENDER_SPLINE_LAST) ? j : j + 1;
            const float* xp = &state.xyz[ord[idx[jp] & 0x3FFFFFFFu]].x;
            const float* xn = &state.xyz[ord[idx[jn] & 0x3FFFFFFFu]].x;
            const float* c4 = &state.xyz[ord[idx[j] & 0x3FFFFFFFu]].x;
            const float* fc = &state.xyz[ord[face[j]]].x;
            float t[3] = { xn[0] - xp[0], xn[1] - xp[1], xn[2] - xp[2] }, F[3] = { fc[0] - c4[0], fc[1] - c4[1], fc[2] - c4[2] };
            float l  = sqrtf(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
            float tl = sqrtf(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]), fl = sqrtf(F[0] * F[0] + F[1] * F[1] + F[2] * F[2]);
            float at = fabsf(a[0] * t[0] + a[1] * t[1] + a[2] * t[2]) / tl, af = fabsf(a[0] * F[0] + a[1] * F[1] + a[2] * F[2]) / fl;
            worst = fmaxf(worst, fmaxf(at, af));
            if (fabsf(l - 1.0f) > 1e-4f || at > 1e-3f || af > 1e-3f) bad++;
        }
        EXPECT_EQ(0, bad);
        printf("    axes: across the chord and the way to the base within %.1e\n", worst);
        free(ax); free(idx); free(face); free(ord);
    }

    /* Picks: of one nucleotide alone (the representation's filter), viewed
       from +z: its ring of six at its centre, the band at its C4' */
    {
        const uint32_t g = nseg / 3;
        const md_nucleic_acid_atoms_t at = NB->segment.atoms[g];
        uint32_t words = (na + 31) / 32;
        uint32_t* bits = calloc(words, 4);
        const md_residue_idx_t res = NB->segment.res_idx[g];
        md_urange_t rr = md_system_residue_atom_range(&sys, (size_t)res);
        for (uint32_t i = (uint32_t)rr.beg; i < (uint32_t)rr.end; ++i) bits[i >> 5] |= 1u << (i & 31);
        rep.visible = f_buf(&fx, bits, (size_t)words * 4);
        ASSERT_EQ(5u, md_render_mol_emit(&topo, dyn, &atoms, &rep, bt, 8));
        /* the ring of six: its atoms, its centre */
        uint32_t* pts = (uint32_t*)read_block(&fx, topo.ring_points[3], (uint64_t)topo.ring_count[3] * 24);
        uint32_t* rs  = (uint32_t*)read_block(&fx, topo.ring_segment[3], (uint64_t)topo.ring_count[3] * 4);
        uint32_t* ord = (uint32_t*)read_block(&fx, topo.atom_order, (uint64_t)topo.point_count * 4);
        int32_t ring = -1;
        for (uint32_t r = 0; r < topo.ring_count[3]; ++r) if (rs[r] == g) ring = (int32_t)r;
        ASSERT_GE(ring, 0);
        float cen[3] = { 0, 0, 0 }, lo[3] = { 1e30f, 1e30f, 1e30f }, hi[3] = { -1e30f, -1e30f, -1e30f };
        for (uint32_t l = 0; l < 6; ++l) for (int k = 0; k < 3; ++k) cen[k] += (&state.xyz[ord[pts[6 * ring + l]]].x)[k] / 6.0f;
        for (uint32_t i = (uint32_t)rr.beg; i < (uint32_t)rr.end; ++i) for (int k = 0; k < 3; ++k) { lo[k] = fminf(lo[k], (&state.xyz[i].x)[k]); hi[k] = fmaxf(hi[k], (&state.xyz[i].x)[k]); }
        const uint32_t W = 128, H = 128;
        float ext = fmaxf(hi[0] - lo[0], hi[1] - lo[1]) * 0.6f + 2.0f;
        cam_t c = { .ortho = true, .half_h = ext, .aspect = 1.0f, .n = 0.0f, .far_ = 1000.0f };
        md_render_camera_t cam = make_camera(&c);
        float ctr[3] = { 0.5f * (lo[0] + hi[0]), 0.5f * (lo[1] + hi[1]), hi[2] + 20.0f };
        cam.view[12] = -ctr[0]; cam.view[13] = -ctr[1]; cam.view[14] = -ctr[2];
        #define PX(x) ((uint32_t)((((x) - ctr[0]) / ext * 0.5f + 0.5f) * W))
        #define PY(y) ((uint32_t)((0.5f - ((y) - ctr[1]) / ext * 0.5f) * H))
        md_render_frame_t f = frame_of(W, H, cam, &bt[3], 1);       /* the rings of six */
        memcpy(f.content, content, sizeof(content));
        ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
        md_render_pick_t pk;
        md_render_mol_hit_t hit;
        ASSERT_TRUE(f_pick(&fx, PX(cen[0]), PY(cen[1]), &pk));
        ASSERT_TRUE(md_render_mol_decode_pick(&pk, &hit));
        EXPECT_EQ(MD_RENDER_MOL_HIT_NUCLEOTIDE, hit.kind);
        EXPECT_EQ(9u, (unsigned)hit.rep_key);
        EXPECT_EQ(g, hit.segment);
        uint32_t near = 0;
        float bd = 1e30f;
        for (uint32_t l = 0; l < 6; ++l) {
            const float* q = &state.xyz[ord[pts[6 * ring + l]]].x;
            float d = (q[0] - hit.position[0]) * (q[0] - hit.position[0]) + (q[1] - hit.position[1]) * (q[1] - hit.position[1]) + (q[2] - hit.position[2]) * (q[2] - hit.position[2]);
            if (d < bd) { bd = d; near = ord[pts[6 * ring + l]]; }
        }
        EXPECT_EQ(near, hit.atom);
        /* the band at the C4': through it */
        f = frame_of(W, H, cam, &bt[1], 1);
        memcpy(f.content, content, sizeof(content));
        f.content[MD_RENDER_CONTENT_GEOMETRY]++;
        ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
        const float* c4 = &state.xyz[at.c4].x;
        ASSERT_TRUE(f_pick(&fx, PX(c4[0]), PY(c4[1]), &pk));
        ASSERT_TRUE(md_render_mol_decode_pick(&pk, &hit));
        EXPECT_EQ(MD_RENDER_MOL_HIT_NUCLEOTIDE, hit.kind);
        EXPECT_EQ(g, hit.segment);
        EXPECT_EQ((uint32_t)at.c4, hit.atom);
        #undef PX
        #undef PY
        free(pts); free(rs); free(ord); free(bits);
    }
    md_gpu_free(fx.gfx, P);
    md_gpu_free(fx.gfx, dyn);
    md_gpu_free(fx.gfx, topo.base);
    md_render_mol_destroy(kernels);
    md_system_free(&sys);
    f_close(&fx);
}

static uint16_t float_to_half(float f) {        /* round to nearest; tiny values to 0 (the tests' are not) */
    uint32_t x;
    memcpy(&x, &f, 4);
    uint32_t sign = (x >> 16) & 0x8000u;
    int32_t  e = (int32_t)((x >> 23) & 0xFFu) - 127 + 15;
    if ((x & 0x7FFFFFFFu) == 0 || e <= 0) return (uint16_t)sign;
    if (e >= 31) return (uint16_t)(sign | 0x7C00u);
    return (uint16_t)(sign | (((uint32_t)e << 10) + (((x & 0x7FFFFFu) + 0x1000u) >> 13)));
}

/* Derived points: the lipids' head centres of the Martini 3 lipid zoo (its
   molecules broken across the box as the run wrote them; the lipids of up
   to 16 beads small molecules, which every regroup moves among the points).
   Each the mean of its head's beads: as given without a cell, at their
   images nearest the first with the frames' cell; still after a regroup;
   its velocity the mean of theirs; standing for its first bead */
/* A derived point as the topology defines it, in double: the weighted sum of
   its members, each at its image nearest the first (L: an orthorhombic cell,
   NULL none) */
typedef struct der_ref_t { uint32_t count; const uint32_t* offset; const uint32_t* member; const float* weight; } der_ref_t;
static void derived_expect(const der_ref_t* D, const float* xyz, const float L[3], uint32_t d, double out[3]) {
    const float* r = xyz + 3 * D->member[D->offset[d]];
    for (int k = 0; k < 3; ++k) out[k] = 0;
    for (uint32_t o = D->offset[d]; o < D->offset[d + 1]; ++o) {
        const float* x = xyz + 3 * D->member[o];
        for (int k = 0; k < 3; ++k) {
            double v = x[k] - r[k];
            if (L) v -= L[k] * round(v / L[k]);
            out[k] += D->weight[o] * (r[k] + v);
        }
    }
}

/* The derived points of the lipids' cartoon (viamd's glyphs) on Martini 3
   lipids: their spheres (a head, a split head as one or its parts), the
   tubes; every derived point where its members' weighted sum is: as given,
   with velocities, after a regroup (the small lipids' beads elsewhere among
   the points), from a frame with the cell and from positions with it (whole
   across the box, each lipid about its first head bead) */
UTEST(render, mol_derived_points) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(load_system(&sys, &state, MD_UNITTEST_DATA_DIR "/martini/martini3_lipids.gro"));
    const uint32_t n = (uint32_t)sys.atom.count;
    const size_t nh = sys.lipid.head.count;
    ASSERT_TRUE(md_unitcell_is_orthorhombic(&state.unitcell));
    const float L[3] = { (float)state.unitcell.x, (float)state.unitcell.y, (float)state.unitcell.z };
    uint64_t content[MD_RENDER_CONTENT_COUNT] = {0};

    md_temp_scope_t temp = md_temp_begin();
    md_render_mol_plan_t plan;
    ASSERT_TRUE(md_render_mol_plan(&plan, &sys, NULL, temp));
    md_render_mol_topo_t topo;
    ASSERT_TRUE(md_render_mol_write(fx.gfx, &topo, &plan, md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, plan.size).gpu, content));
    const uint32_t nd = topo.derived_count;
    der_ref_t D = { nd, NULL, NULL, NULL };
    {
        const uint8_t* img = plan.data;
        const uint32_t* off = (const uint32_t*)(img + plan.layout.derived_offset);
        uint32_t* o2 = malloc(((size_t)nd + 1) * 4); memcpy(o2, off, ((size_t)nd + 1) * 4);
        uint32_t* m2 = malloc((size_t)off[nd] * 4 + 4); memcpy(m2, img + plan.layout.derived_member, (size_t)off[nd] * 4);
        float*    w2 = malloc((size_t)off[nd] * 4 + 4); memcpy(w2, img + plan.layout.derived_weight, (size_t)off[nd] * 4);
        D.offset = o2; D.member = m2; D.weight = w2;
    }
    md_temp_end(temp);
    uint32_t small_heads = 0, split = 0;
    for (size_t h = 0; h < nh; ++h) {
        const md_atom_idx_t* at = md_lipid_head_atoms(&sys.lipid, h);
        int32_t m = -1;
        for (size_t k = 0; k < sys.molecule.count && m < 0; ++k) {
            for (uint32_t j = sys.molecule.offset[k]; j < sys.molecule.offset[k + 1]; ++j) if (sys.molecule.atom_idx[j] == at[0]) m = (int32_t)k;
        }
        if (m >= 0 && sys.molecule.offset[m + 1] - sys.molecule.offset[m] <= MD_RENDER_MOL_SMALL) small_heads++;
    }
    for (uint32_t d = 0; d < nd; ++d) {
        double a[3], b[3];
        derived_expect(&D, &state.xyz[0].x, NULL, d, a);
        derived_expect(&D, &state.xyz[0].x, L, d, b);
        split += fabs(a[0] - b[0]) + fabs(a[1] - b[1]) + fabs(a[2] - b[2]) > 1.0;
    }
    const uint32_t* sc = topo.lipid_sphere_count;
    printf("    %zu lipids, %zu heads (%u of small molecules): %u derived points after %u atoms' points, %u split by the box; spheres %u + %u (or %u parts), tubes' control points %u %u %u\n",
           sys.lipid.count, nh, small_heads, nd, topo.point_count, split, sc[0], sc[1], sc[2], topo.lipid_tube_count[0], topo.lipid_tube_count[1], topo.lipid_tube_count[2]);
    EXPECT_EQ((uint32_t)nh, sc[0] + sc[1]);        /* a sphere per head (every lipid's heads known here) */
    EXPECT_GT(sc[2], sc[1]);                        /* the split heads' parts */
    EXPECT_GT(topo.lipid_tube_count[0], 100u);
    EXPECT_GT(topo.lipid_tube_count[2], 0u);        /* a cholesterol's slab */
    EXPECT_GT(small_heads, 0u);
    EXPECT_GT(split, 0u);
    {   /* a head's sphere: its radius viamd's, standing for its heaviest (a bead: the phosphate's, or its first) */
        uint32_t* sph = (uint32_t*)read_block(&fx, topo.lipid_spheres, ((size_t)sc[0] + sc[1] + sc[2]) * 4);
        float* rad = (float*)read_block(&fx, topo.base + plan.layout.derived_radius, (size_t)nd * 4);
        int bad = 0;
        for (uint32_t k = 0; k < sc[0] + sc[1] + sc[2]; ++k) {
            const uint32_t d = sph[k] - topo.point_count;
            const uint32_t m = D.offset[d + 1] - D.offset[d] - (D.weight[D.offset[d]] == 0.0f ? 1u : 0u);   /* less a lead of weight 0 */
            float r = 1.6f * cbrtf(4.0f * (float)m);
            r = r < 1.8f ? 1.8f : r > 7.0f ? 7.0f : r;
            bad += d >= nd || fabsf(rad[d] - r) > 1e-5f;
        }
        EXPECT_EQ(0, bad);
        free(sph); free(rad);
    }

    md_render_mol_ctx_t kernels = md_render_mol_create(fx.dev, NULL);
    ASSERT_TRUE(kernels != NULL);
    md_gpu_addr_t dyn = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, topo.dynamic_size_velocity).gpu;
    const uint32_t np = topo.point_count;
    #define CHECK_DERIVED(cellL, tol, label) do { \
        float* pts = (float*)read_block(&fx, dyn, ((size_t)np + nd) * 16); \
        double worst = 0; \
        for (uint32_t d = 0; d < nd; ++d) { \
            double e[3]; derived_expect(&D, &state.xyz[0].x, cellL, d, e); \
            for (int k = 0; k < 3; ++k) worst = fmax(worst, fabs(pts[4 * (np + d) + k] - e[k])); \
        } \
        printf("    %s: the derived points within %.2e A of their sums\n", label, worst); \
        EXPECT_LT(worst, tol); \
        free(pts); \
    } while (0)

    /* as given, no cell; with velocities (packed half3) */
    md_gpu_addr_t P = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, (uint64_t)n * 12).gpu;
    ASSERT_TRUE(md_render_mol_upload_state(fx.gfx, P, &state, content));
    uint16_t* vel = malloc((size_t)n * 8);
    float* velf = malloc((size_t)n * 12);
    for (uint32_t i = 0; i < n; ++i) {
        for (int d = 0; d < 3; ++d) { velf[3 * i + d] = (float)((int)((i * 7u + (uint32_t)d * 13u) % 29u) - 14) * 0.25f; vel[4 * i + d] = float_to_half(velf[3 * i + d]); }
        vel[4 * i + 3] = 0;
    }
    md_render_mol_atoms_t atoms = {0};
    atoms.position = P;
    atoms.velocity = f_buf(&fx, vel, (size_t)n * 8);
    ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, 0, content));
    CHECK_DERIVED(NULL, 1e-4, "as given");
    {
        uint16_t* pv = (uint16_t*)read_block(&fx, dyn + (topo.dynamic_size), ((size_t)np + nd) * 8);
        double worst = 0;
        for (uint32_t d = 0; d < nd; ++d) {
            for (int k = 0; k < 3; ++k) {
                double e = 0;
                for (uint32_t o = D.offset[d]; o < D.offset[d + 1]; ++o) e += D.weight[o] * velf[3 * D.member[o] + k];
                worst = fmax(worst, fabs(half_to_float(pv[4 * (np + d) + k]) - e));
            }
        }
        printf("    velocities: the derived points within %.2e of their sums\n", worst);
        EXPECT_LT(worst, 1e-2);
        free(pv);
    }
    /* regrouped: the small lipids' beads elsewhere among the points, the derived where they were */
    uint32_t* ord0 = (uint32_t*)read_block(&fx, topo.atom_order, (size_t)np * 4);
    ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, MD_RENDER_MOL_UPDATE_REGROUP, content));
    uint32_t* ord1 = (uint32_t*)read_block(&fx, topo.atom_order, (size_t)np * 4);
    EXPECT_NE(0, memcmp(ord0, ord1, (size_t)np * 4));
    CHECK_DERIVED(NULL, 1e-4, "after a regroup");
    free(ord0); free(ord1);

    /* from positions with their cell: whole */
    {
        md_render_mol_atoms_t ac = atoms;
        md_render_mol_frame_cell(ac.cell, &state.unitcell);
        ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &ac, 0, content));
        CHECK_DERIVED(L, 1e-4, "with the positions' cell, whole");
    }
    /* from a frame with the cell: whole */
    {
        float cell[9];
        md_render_mol_frame_cell(cell, &state.unitcell);
        uint16_t* q = malloc((size_t)n * 6);
        md_render_mol_frame_t f;
        md_render_mol_frame_quantize(&f, q, &state.xyz[0].x, 0, n, cell, NULL);
        md_render_mol_frames_t fr;
        ASSERT_TRUE(md_render_mol_frames_init(&fr, &topo, md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, md_render_mol_frames_size(&topo, 1)).gpu, 1));
        ASSERT_TRUE(md_render_mol_frame_upload(kernels, fx.gfx, &topo, dyn, &fr, 0, &f, MD_RENDER_MOL_UPDATE_REGROUP));
        md_render_mol_atoms_t af = {0};
        af.frames = &fr;
        af.frame_count = 1;
        af.frame_w[0] = 1.0f;
        ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &af, 0, content));
        const double qtol = 2.0 * fmax(fmax(f.scale[0], f.scale[1]), f.scale[2]);
        CHECK_DERIVED(L, qtol, "from a frame, whole");
        free(q);
    }
    #undef CHECK_DERIVED
    free((void*)D.offset); free((void*)D.member); free((void*)D.weight);
    free(vel); free(velf);
    md_render_mol_destroy(kernels);
    md_system_free(&sys);
    f_close(&fx);
}

/* The bonds that are not covalent dashed, in their kind's style (md_gl's
   BOND DASHES, the default table): the coordinate bonds of a catalyst. Their
   batches' styles by representation and by the caller's table (radius,
   dashes, colour; radius 0 not drawn; a negative length refused); drawn,
   the dashes cover part of what the solid sticks do, never more */
UTEST(render, mol_bonds_dashed_by_kind) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(load_system(&sys, &state, MD_UNITTEST_DATA_DIR "/catalyst.gro"));
    uint64_t content[MD_RENDER_CONTENT_COUNT] = {0};
    md_temp_scope_t temp = md_temp_begin();
    md_render_mol_plan_t plan;
    ASSERT_TRUE(md_render_mol_plan(&plan, &sys, NULL, temp));
    md_render_mol_topo_t topo;
    ASSERT_TRUE(md_render_mol_write(fx.gfx, &topo, &plan, md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, plan.size).gpu, content));
    md_temp_end(temp);
    const uint32_t ncoord = topo.bond_kind_offset[MD_RENDER_MOL_BOND_KIND(MD_BOND_KIND_COORDINATE) + 1] - topo.bond_kind_offset[MD_RENDER_MOL_BOND_KIND(MD_BOND_KIND_COORDINATE)];
    ASSERT_GT(ncoord, 0u);
    md_render_mol_ctx_t kernels = md_render_mol_create(fx.dev, NULL);
    ASSERT_TRUE(kernels != NULL);
    md_gpu_addr_t dyn = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, topo.dynamic_size).gpu;
    md_gpu_addr_t P = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, (uint64_t)sys.atom.count * 12).gpu;
    ASSERT_TRUE(md_render_mol_upload_state(fx.gfx, P, &state, content));
    md_render_mol_atoms_t atoms = {0};
    atoms.position = P;
    ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, 0, content));

    const float r = 0.2f;
    md_render_mol_rep_t rep = {0};
    rep.key = 5; rep.type = MD_RENDER_MOL_REP_LICORICE; rep.params.licorice.radius = r;
    md_render_batch_t bt[8];
    ASSERT_EQ(3u, md_render_mol_emit(&topo, dyn, &atoms, &rep, bt, 8));     /* covalent, coordinate, atoms */
    const md_render_capsule_style_t* cv = &bt[0].style.capsules;
    const md_render_capsule_style_t* cd = &bt[1].style.capsules;
    EXPECT_EQ((uint64_t)(5 << 8 | (MD_RENDER_MOL_PART_BONDS + MD_RENDER_MOL_BOND_KIND(MD_BOND_KIND_COORDINATE))), bt[1].key);
    EXPECT_EQ(0.0f, cv->dash_length);
    EXPECT_NEAR(0.70f * r, cd->radius, 1e-6f);
    EXPECT_NEAR(2.0f * r, cd->dash_length, 1e-6f);
    EXPECT_NEAR(0.8f * r, cd->gap_length, 1e-6f);
    EXPECT_EQ(0.0f, cd->dash_ends.radius_scale);       /* licorice: its joints, r */
    EXPECT_NEAR(r, cd->dash_ends.radius_add, 1e-6f);
    EXPECT_EQ((int)MD_RENDER_COLOR_POINTS, (int)bt[1].color.mode);   /* the atoms' colours */
    md_render_batch_t dashed = bt[1];
    const uint32_t kc = MD_RENDER_MOL_BOND_KIND(MD_BOND_KIND_COORDINATE), kv = MD_RENDER_MOL_BOND_KIND(MD_BOND_KIND_COVALENT);
    md_render_mol_bond_style_t tab[MD_RENDER_MOL_BOND_KIND_COUNT];
    memcpy(tab, md_render_mol_bond_styles, sizeof(tab));
    tab[kc] = (md_render_mol_bond_style_t){ 1.0f, 0.0f, 0.0f, 0xFF00FF00u };     /* solid, a colour of its own */
    rep.bond_style = tab;
    ASSERT_EQ(3u, md_render_mol_emit(&topo, dyn, &atoms, &rep, bt, 8));
    EXPECT_EQ(0.0f, bt[1].style.capsules.dash_length);
    EXPECT_NEAR(r, bt[1].style.capsules.radius, 1e-6f);
    EXPECT_EQ((int)MD_RENDER_COLOR_UNIFORM, (int)bt[1].color.mode);
    EXPECT_EQ(0xFF00FF00u, bt[1].color.uniform);
    EXPECT_EQ((int)MD_RENDER_COLOR_POINTS, (int)bt[0].color.mode);
    md_render_batch_t solid = bt[1];
    solid.color = bt[0].color;                                                 /* compared below by coverage only */
    tab[kv] = (md_render_mol_bond_style_t){ 1.0f, 2.4f, 0.6f, 0 };             /* md_gl's DASHES_ALL for the covalent */
    ASSERT_EQ(3u, md_render_mol_emit(&topo, dyn, &atoms, &rep, bt, 8));
    EXPECT_NEAR(2.4f * r, bt[0].style.capsules.dash_length, 1e-6f);
    EXPECT_NEAR(0.6f * r, bt[0].style.capsules.gap_length, 1e-6f);
    EXPECT_NEAR(1.0f * r, bt[0].style.capsules.radius, 1e-6f);
    tab[kc].radius = 0.0f;                                                     /* not drawn */
    ASSERT_EQ(2u, md_render_mol_emit(&topo, dyn, &atoms, &rep, bt, 8));
    EXPECT_EQ((uint64_t)(5 << 8 | (MD_RENDER_MOL_PART_BONDS + kv)), bt[0].key);
    EXPECT_EQ((uint64_t)(5 << 8 | MD_RENDER_MOL_PART_ATOMS), bt[1].key);
    tab[kc] = (md_render_mol_bond_style_t){ 0.7f, -1.0f, 0.8f, 0 };
    EXPECT_EQ(0u, md_render_mol_emit(&topo, dyn, &atoms, &rep, bt, 8));      /* refused */
    rep.bond_style = NULL;
    md_render_mol_rep_t bs = {0};
    bs.key = 6; bs.type = MD_RENDER_MOL_REP_BALL_AND_STICK; bs.params.ball_and_stick.ball_scale = 0.3f; bs.params.ball_and_stick.stick_radius = 0.12f;
    ASSERT_EQ(3u, md_render_mol_emit(&topo, dyn, &atoms, &bs, bt, 8));
    EXPECT_NEAR(0.3f, bt[1].style.capsules.dash_ends.radius_scale, 1e-6f);   /* its balls */
    EXPECT_EQ(0.0f, bt[1].style.capsules.dash_ends.radius_add);
    EXPECT_NEAR(2.0f * 0.12f, bt[1].style.capsules.dash_length, 1e-6f);

    /* drawn alone: the dashes within the solid sticks, and less of them */
    float lo[3] = { 1e30f, 1e30f, 1e30f }, hi[3] = { -1e30f, -1e30f, -1e30f };
    for (size_t i = 0; i < sys.atom.count; ++i) for (int k = 0; k < 3; ++k) { lo[k] = fminf(lo[k], (&state.xyz[i].x)[k]); hi[k] = fmaxf(hi[k], (&state.xyz[i].x)[k]); }
    const uint32_t W = 256, H = 256;
    float ext = 0.55f * fmaxf(hi[0] - lo[0], hi[1] - lo[1]) + 1.0f;
    cam_t c = { .ortho = true, .half_h = ext, .aspect = 1.0f, .n = 0.0f, .far_ = 1000.0f };
    md_render_camera_t cam = make_camera(&c);
    cam.view[12] = -0.5f * (lo[0] + hi[0]); cam.view[13] = -0.5f * (lo[1] + hi[1]); cam.view[14] = -(hi[2] + 20.0f);
    uint32_t* ids[2];
    md_render_targets_t tg = md_render_view_targets(fx.view);
    const md_render_batch_t* which[2] = { &solid, &dashed };
    for (int k = 0; k < 2; ++k) {
        md_render_frame_t f = frame_of(W, H, cam, which[k], 1);
        memcpy(f.content, content, sizeof(content));
        f.content[MD_RENDER_CONTENT_GEOMETRY] += (uint64_t)k + 1;
        ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
        ids[k] = malloc(W * H * 8);
        tg = md_render_view_targets(fx.view);       /* of this size */
        ASSERT_TRUE(f_read(&fx, tg.id, ids[k], W * H * 8));
    }
    int ns = 0, nd = 0, outside = 0;
    for (uint32_t i = 0; i < W * H; ++i) {
        ns += ids[0][2 * i + 1] != 0;
        if (ids[1][2 * i + 1] == 0) continue;
        nd++;
        outside += ids[0][2 * i + 1] == 0;      /* (where sticks cross, another may be in front) */
    }
    printf("    %u coordinate bonds: %d pixels solid, %d dashed\n", ncoord, ns, nd);
    EXPECT_EQ(0, outside);
    EXPECT_GT(nd, ns / 5);
    EXPECT_LT(nd, ns * 3 / 4);
    free(ids[0]); free(ids[1]);
    md_render_mol_destroy(kernels);
    md_system_free(&sys);
    f_close(&fx);
}

/* md_gl's chain gaps on 9gfu (residues missing in its chains, as md_mmcif's
   test has it): the cartoon's from CA to CA, the links' from the C to the N;
   emitted by the representations that ask for them; along the cartoon's,
   dashes and gaps, a dash picking as the residue before the break and the
   end nearer the hit */
UTEST(render, mol_chain_gaps) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_mmcif_system_init_from_file(&sys, &state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/9gfu.cif"), NULL));
    ASSERT_TRUE(md_system_infer(&sys, &state, NULL));
    uint32_t brk[64], nbrk = 0;
    for (uint32_t i = 0; i + 1 < sys.residue.count && nbrk < 64; ++i) {
        if (md_system_residue_flags(&sys, i) & MD_RESIDUE_FLAG_BREAK_AFTER) brk[nbrk++] = i;
    }
    ASSERT_GE(nbrk, 4u);
    /* each residue's CA */
    const md_protein_backbone_data_t* PB = &sys.protein_backbone;
    int32_t* ca = malloc(sys.residue.count * sizeof(int32_t));
    for (size_t i = 0; i < sys.residue.count; ++i) ca[i] = -1;
    for (size_t g = 0; g < PB->segment.count; ++g) if (PB->segment.res_idx[g] >= 0) ca[PB->segment.res_idx[g]] = PB->segment.atoms[g].ca;

    md_temp_scope_t temp = md_temp_begin();
    md_render_mol_plan_t plan;
    ASSERT_TRUE(md_render_mol_plan(&plan, &sys, NULL, temp));
    const md_render_mol_topo_t* L = &plan.layout;
    EXPECT_EQ(nbrk, L->gap_count[0]);
    EXPECT_EQ(nbrk, L->gap_count[1]);
    const uint8_t*  img   = plan.data;
    const uint32_t* order = (const uint32_t*)(img + L->atom_order);
    const uint32_t* gp    = (const uint32_t*)(img + L->gap_points);
    const uint32_t* gr    = (const uint32_t*)(img + L->gap_residue);
    int wrong = 0;
    for (uint32_t k = 0; k < nbrk && k < L->gap_count[0] && k < L->gap_count[1]; ++k) {
        const uint32_t i = brk[k];
        md_atom_idx_t la[2];
        ASSERT_TRUE(md_system_chain_break_atoms(la, &sys, i));
        wrong += gr[k] != i || (int32_t)order[gp[2 * k]] != ca[i] || (int32_t)order[gp[2 * k + 1]] != ca[i + 1];
        const uint32_t j = L->gap_count[0] + k;
        wrong += gr[j] != i || (int32_t)order[gp[2 * j]] != la[0] || (int32_t)order[gp[2 * j + 1]] != la[1];
    }
    EXPECT_EQ(0, wrong);

    /* A nucleic acid's from C4' to C4': 8g7u without the sixth nucleotide of
       its RNA (chain E) */
    {
        FILE* fp = fopen(MD_UNITTEST_DATA_DIR "/8g7u.cif", "rb");
        ASSERT_TRUE(fp != NULL);
        fseek(fp, 0, SEEK_END);
        const long fsz = ftell(fp);
        fseek(fp, 0, SEEK_SET);
        char* txt = malloc((size_t)fsz + 1), *cut = malloc((size_t)fsz + 1);
        ASSERT_EQ((size_t)fsz, fread(txt, 1, (size_t)fsz, fp));
        fclose(fp);
        txt[fsz] = 0;
        size_t m = 0, dropped = 0;
        for (char* ln = txt; *ln; ) {
            char* nl = strchr(ln, '\n');
            const size_t ll = nl ? (size_t)(nl - ln) + 1 : strlen(ln);
            char f[10][16] = {{0}};
            if (!strncmp(ln, "ATOM", 4)) sscanf(ln, "%15s %15s %15s %15s %15s %15s %15s %15s %15s", f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8]);
            if (!strcmp(f[6], "E") && !strcmp(f[8], "6")) dropped++;
            else { memcpy(cut + m, ln, ll); m += ll; }
            ln += ll;
        }
        cut[m] = 0;
        EXPECT_GT(dropped, 10u);
        md_allocator_i* a2 = md_vm_arena_create(MD_GIGABYTES(1));
        md_system_t ns = { .alloc = a2 };
        md_system_state_t nst = { .alloc = a2 };
        ASSERT_TRUE(md_mmcif_system_init_from_str(&ns, &nst, (md_str_t){ cut, m }, NULL));
        ASSERT_TRUE(md_system_infer(&ns, &nst, NULL));
        int32_t* c4 = malloc((ns.residue.count + 1) * sizeof(int32_t));
        for (size_t i = 0; i < ns.residue.count; ++i) c4[i] = -1;
        const md_nucleic_backbone_data_t* NB = &ns.nucleic_backbone;
        for (size_t g = 0; g < NB->segment.count; ++g) if (NB->segment.res_idx[g] >= 0) c4[NB->segment.res_idx[g]] = NB->segment.atoms[g].c4;
        md_temp_scope_t t2 = md_temp_begin();
        md_render_mol_plan_t np;
        ASSERT_TRUE(md_render_mol_plan(&np, &ns, NULL, t2));
        const md_render_mol_topo_t* NL = &np.layout;
        const uint32_t* no = (const uint32_t*)((const uint8_t*)np.data + NL->atom_order);
        const uint32_t* ngp = (const uint32_t*)((const uint8_t*)np.data + NL->gap_points);
        const uint32_t* ngr = (const uint32_t*)((const uint8_t*)np.data + NL->gap_residue);
        int nucleic = 0;
        for (uint32_t k = 0; k < NL->gap_count[0]; ++k) {
            const uint32_t i = ngr[k];
            if (c4[i] < 0 && c4[i + 1] < 0) continue;
            nucleic++;
            EXPECT_EQ(c4[i], (int32_t)no[ngp[2 * k]]);
            EXPECT_EQ(c4[i + 1], (int32_t)no[ngp[2 * k + 1]]);
            md_atom_idx_t la[2];
            EXPECT_TRUE(md_system_chain_break_atoms(la, &ns, i));
            bool link = false;          /* and its link, from the O3' to the P */
            for (uint32_t j = NL->gap_count[0]; j < NL->gap_count[0] + NL->gap_count[1]; ++j) {
                link |= ngr[j] == i && (int32_t)no[ngp[2 * j]] == la[0] && (int32_t)no[ngp[2 * j + 1]] == la[1];
            }
            EXPECT_TRUE(link);
            EXPECT_TRUE(md_str_eq_cstr(md_atom_name(&ns.atom, (size_t)la[0]), "O3'"));
            EXPECT_TRUE(md_str_eq_cstr(md_atom_name(&ns.atom, (size_t)la[1]), "P"));
        }
        printf("    8g7u without a nucleotide: %u gaps of its cartoon, %d of them a nucleic acid's\n", NL->gap_count[0], nucleic);
        EXPECT_EQ(1, nucleic);
        md_temp_end(t2);
        free(c4); free(txt); free(cut);
        md_system_free(&ns);
        md_vm_arena_destroy(a2);
    }

    fix_t fx;
    if (!f_open(&fx)) { md_temp_end(temp); free(ca); md_system_free(&sys); UTEST_SKIP(f_skip_reason()); }
    uint64_t content[MD_RENDER_CONTENT_COUNT] = {0};
    md_render_mol_topo_t topo;
    ASSERT_TRUE(md_render_mol_write(fx.gfx, &topo, &plan, md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, plan.size).gpu, content));
    md_temp_end(temp);
    md_render_mol_ctx_t kernels = md_render_mol_create(fx.dev, NULL);
    ASSERT_TRUE(kernels != NULL);
    md_gpu_addr_t dyn = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, topo.dynamic_size).gpu;
    md_gpu_addr_t P = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, (uint64_t)sys.atom.count * 12).gpu;
    ASSERT_TRUE(md_render_mol_upload_state(fx.gfx, P, &state, content));
    md_render_mol_atoms_t atoms = {0};
    atoms.position = P;
    ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, 0, content));

    /* who emits them */
    md_render_mol_rep_t rep = {0};
    rep.key = 9; rep.type = MD_RENDER_MOL_REP_CARTOON;
    md_render_batch_t bt[16];
    const uint32_t n0 = md_render_mol_emit(&topo, dyn, &atoms, &rep, bt, 16);
    rep.chain_gaps = true;
    const uint32_t n1 = md_render_mol_emit(&topo, dyn, &atoms, &rep, bt, 16);
    ASSERT_EQ(n0 + 1, n1);
    md_render_batch_t gaps = bt[n1 - 1];
    EXPECT_EQ((uint64_t)(9 << 8 | MD_RENDER_MOL_PART_GAPS), gaps.key);
    EXPECT_EQ(nbrk, gaps.count);
    EXPECT_NEAR(0.14f, gaps.style.capsules.radius, 1e-6f);
    EXPECT_NEAR(0.56f, gaps.style.capsules.dash_length, 1e-6f);
    EXPECT_NEAR(0.56f, gaps.style.capsules.gap_length, 1e-6f);
    EXPECT_EQ(0.0f, gaps.style.capsules.dash_ends.radius_add);
    EXPECT_EQ(0.0f, gaps.color.sharpness);
    md_render_mol_rep_t lic = {0};
    lic.key = 10; lic.type = MD_RENDER_MOL_REP_LICORICE; lic.params.licorice.radius = 0.2f;
    const uint32_t m0 = md_render_mol_emit(&topo, dyn, &atoms, &lic, bt, 16);
    lic.chain_gaps = true;
    ASSERT_EQ(m0 + 1, md_render_mol_emit(&topo, dyn, &atoms, &lic, bt, 16));
    EXPECT_EQ((uint64_t)(10 << 8 | MD_RENDER_MOL_PART_GAPS), bt[m0 - 1].key);       /* before the atoms */
    EXPECT_EQ(topo.gap_points + (md_gpu_addr_t)topo.gap_count[0] * 8u, bt[m0 - 1].index);
    EXPECT_NEAR(0.1f, bt[m0 - 1].style.capsules.radius, 1e-6f);
    EXPECT_NEAR(0.4f, bt[m0 - 1].style.capsules.dash_length, 1e-6f);
    EXPECT_NEAR(0.2f, bt[m0 - 1].style.capsules.dash_ends.radius_add, 1e-6f);   /* from the joints' spheres */
    md_render_mol_rep_t sf = {0};
    sf.key = 11; sf.type = MD_RENDER_MOL_REP_SPACEFILL; sf.params.spacefill.radius_scale = 1.0f; sf.chain_gaps = true;
    EXPECT_EQ(1u, md_render_mol_emit(&topo, dyn, &atoms, &sf, bt, 16));          /* no sticks, no gaps */

    /* The first gap side on, along x: picks along it */
    const float* pa = &state.xyz[ca[brk[0]]].x;
    const float* pb = &state.xyz[ca[brk[0] + 1]].x;
    double x[3] = { pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2] }, len = sqrt(dot3(x, x));
    for (int k = 0; k < 3; ++k) x[k] /= len;
    double z[3] = { 0, 0, 1 };
    if (fabs(x[2]) > 0.9) { z[0] = 1; z[2] = 0; }
    double y[3];
    cross3(z, x, y);
    double ly = sqrt(dot3(y, y));
    for (int k = 0; k < 3; ++k) y[k] /= ly;
    cross3(x, y, z);
    double eye[3];
    for (int k = 0; k < 3; ++k) eye[k] = 0.5 * (pa[k] + pb[k]) + 50.0 * z[k];
    const uint32_t W = 256, H = 64;
    cam_t c = { .ortho = true, .half_h = (float)(0.55 * len / 4.0), .aspect = 4.0f, .n = 0.0f, .far_ = 1000.0f };
    md_render_camera_t cam = make_camera(&c);
    for (int k = 0; k < 3; ++k) {
        cam.view[4 * k + 0] = (float)x[k];
        cam.view[4 * k + 1] = (float)y[k];
        cam.view[4 * k + 2] = (float)z[k];
    }
    cam.view[12] = (float)-dot3(x, eye); cam.view[13] = (float)-dot3(y, eye); cam.view[14] = (float)-dot3(z, eye);
    md_render_frame_t f = frame_of(W, H, cam, &gaps, 1);
    memcpy(f.content, content, sizeof(content));
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    int hits = 0, misses = 0, bad = 0, runs = 0;
    bool was = false;
    for (int s = 3; s <= 97; ++s) {
        const double t = s / 100.0, sx = (t - 0.5) * len;
        const uint32_t px = (uint32_t)((sx / (c.half_h * c.aspect) * 0.5 + 0.5) * W);
        md_render_pick_t pk;
        md_render_mol_hit_t hit;
        bool on = f_pick(&fx, px, H / 2, &pk) && md_render_mol_decode_pick(&pk, &hit) && hit.kind == MD_RENDER_MOL_HIT_CHAIN_GAP && hit.residue == brk[0];
        if (on) {
            hits++;
            if (fabs(t - 0.5) > 0.02) bad += hit.atom != (uint32_t)(t < 0.5 ? ca[brk[0]] : ca[brk[0] + 1]);
            bad += hit.rep_key != 9;
        } else {
            misses++;
        }
        runs += on && !was;
        was = on;
    }
    const double D = 0.56, nd = fmax(2.0, floor((len + D) / (2 * D) + 0.5));
    printf("    %u breaks; the first %.2f A: %d of %d picks on a dash, %d dashes met (%.0f drawn)\n", nbrk, len, hits, hits + misses, runs, nd);
    EXPECT_EQ(0, bad);
    EXPECT_GT(hits, 30);
    EXPECT_GT(misses, 30);
    EXPECT_GE(runs, (int)nd - 1);
    md_render_mol_destroy(kernels);
    free(ca);
    md_system_free(&sys);
    f_close(&fx);
}

/* Oct-7-ynal, built by hand (ideal lengths and angles, no file of it): a
   triple bond, a double, 21 atoms (not a small molecule) */
static const char octynal_pdb[] =
    "HETATM    1 H1   OYN A   1      -1.060   0.000   0.000  1.00  0.00           H\n"
    "HETATM    2 C2   OYN A   1       0.000   0.000   0.000  1.00  0.00           C\n"
    "HETATM    3 C3   OYN A   1       1.200   0.000   0.000  1.00  0.00           C\n"
    "HETATM    4 C4   OYN A   1       2.660   0.000   0.000  1.00  0.00           C\n"
    "HETATM    5 C5   OYN A   1       3.921   0.867   0.000  1.00  0.00           C\n"
    "HETATM    6 C6   OYN A   1       5.182   0.000   0.000  1.00  0.00           C\n"
    "HETATM    7 C7   OYN A   1       6.443   0.867   0.000  1.00  0.00           C\n"
    "HETATM    8 C8   OYN A   1       7.704   0.000   0.000  1.00  0.00           C\n"
    "HETATM    9 C9   OYN A   1       8.948   0.855   0.000  1.00  0.00           C\n"
    "HETATM   10 H10  OYN A   1       2.847  -0.601  -0.890  1.00  0.00           H\n"
    "HETATM   11 H11  OYN A   1       2.847  -0.601   0.890  1.00  0.00           H\n"
    "HETATM   12 H12  OYN A   1       3.921   1.496   0.890  1.00  0.00           H\n"
    "HETATM   13 H13  OYN A   1       3.921   1.496  -0.890  1.00  0.00           H\n"
    "HETATM   14 H14  OYN A   1       5.182  -0.629  -0.890  1.00  0.00           H\n"
    "HETATM   15 H15  OYN A   1       5.182  -0.629   0.890  1.00  0.00           H\n"
    "HETATM   16 H16  OYN A   1       6.443   1.496   0.890  1.00  0.00           H\n"
    "HETATM   17 H17  OYN A   1       6.443   1.496  -0.890  1.00  0.00           H\n"
    "HETATM   18 H18  OYN A   1       7.704  -0.629  -0.890  1.00  0.00           H\n"
    "HETATM   19 H19  OYN A   1       7.704  -0.629   0.890  1.00  0.00           H\n"
    "HETATM   20 O20  OYN A   1      10.040   0.334   0.000  1.00  0.00           O\n"
    "HETATM   21 H21  OYN A   1       8.862   1.952   0.000  1.00  0.00           H\n"
    "CONECT    1    2\n"
    "CONECT    2    1    3\n"
    "CONECT    3    2    4\n"
    "CONECT    4    3    5   10   11\n"
    "CONECT    5    4    6   12   13\n"
    "CONECT    6    5    7   14   15\n"
    "CONECT    7    6    8   16   17\n"
    "CONECT    8    7    9   18   19\n"
    "CONECT    9    8   20   21\n"
    "CONECT   10    4\n"
    "CONECT   11    4\n"
    "CONECT   12    5\n"
    "CONECT   13    5\n"
    "CONECT   14    6\n"
    "CONECT   15    6\n"
    "CONECT   16    7\n"
    "CONECT   17    7\n"
    "CONECT   18    8\n"
    "CONECT   19    8\n"
    "CONECT   20    9\n"
    "CONECT   21    9\n"
    "END\n";

/* The double bond's side as md_render_mol's bond_sides has it, in double */
static void ref_end_lines(const double ax[3], const double p[3], const double* c0, const double* c1, const double u[3], const double v[3], double s[2]) {
    const double* cs[2] = { c0, c1 };
    double acc[2] = { 0, 0 };
    for (int k = 0; k < 2; ++k) {
        if (!cs[k]) continue;
        double c[3] = { cs[k][0] - p[0], cs[k][1] - p[1], cs[k][2] - p[2] }, cc = dot3(c, c);
        if (!(cc > 0 && cc < 16.0)) { cs[k] = NULL; continue; }
        double x = dot3(c, u), y = dot3(c, v);
        acc[0] += (x * x - y * y) / cc; acc[1] += 2 * x * y / cc;
    }
    if (cs[0] && cs[1]) {
        double a[3] = { cs[0][0] - p[0], cs[0][1] - p[1], cs[0][2] - p[2] }, b[3] = { cs[1][0] - p[0], cs[1][1] - p[1], cs[1][2] - p[2] }, m[3];
        cross3(a, b, m);
        double so = fabs(dot3(ax, m)) / sqrt(fmax(dot3(m, m), 1e-20));
        double t = so <= 0.17 ? 0 : so >= 0.5 ? 1 : (so - 0.17) / (0.5 - 0.17);
        t = t * t * (3 - 2 * t);
        acc[0] *= 1 - t; acc[1] *= 1 - t;
    }
    s[0] += acc[0]; s[1] += acc[1];
}

/* md_gl's bond orders on 1LAF: the double bonds of the plan are those of
   order 2, neither aromatic nor delocalized, the sticks of the covalent
   class leave them out; the partners those of every delocalized bond, small
   molecules' (ions) too; the update's sides are the reference's; licorice
   with bond_orders emits the two thinner sticks, picked as the bond */
UTEST(render, mol_bond_orders) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(load_system(&sys, &state, MD_UNITTEST_DATA_DIR "/1LAF.pdb"));
    md_temp_scope_t temp = md_temp_begin();
    md_render_mol_plan_t plan;
    ASSERT_TRUE(md_render_mol_plan(&plan, &sys, NULL, temp));
    const md_render_mol_topo_t* L = &plan.layout;
    const uint8_t*  img   = plan.data;
    const uint32_t* order = (const uint32_t*)(img + L->atom_order);
    const uint32_t* op    = (const uint32_t*)(img + L->orders.order_atoms);
    const uint32_t* ob    = (const uint32_t*)(img + L->orders.order_id);
    const uint32_t* orf   = (const uint32_t*)(img + L->orders.order_ref);
    const uint32_t  nd = L->orders.order_count[0], nt = L->orders.order_count[1];
    ASSERT_GT(nd, 100u);
    int wrong = 0, both_ends = 0;
    uint8_t* is_ord = calloc(sys.bond.count + 1, 1);
    for (uint32_t k = 0; k < nd + nt; ++k) {
        const uint32_t bi = ob[k], f = sys.bond.flags[bi];
        is_ord[bi] = 1;
        wrong += md_bond_order(f) != (k < nd ? 2 : 3) || (f & (MD_BOND_FLAG_AROMATIC | MD_BOND_FLAG_DELOCALIZED)) || md_bond_kind(f) != MD_BOND_KIND_COVALENT;
        const uint32_t a = op[2 * k], b = op[2 * k + 1];
        wrong += !((a == (uint32_t)sys.bond.pairs[bi].idx[0] && b == (uint32_t)sys.bond.pairs[bi].idx[1]));
        if (k >= nd) continue;
        int ends = 0;
        for (int e = 0; e < 2; ++e) {
            const uint32_t at = e ? b : a, other = e ? a : b;
            bool any = false;
            for (int q = 0; q < 2; ++q) {
                const uint32_t r = orf[4 * k + 2 * e + q];
                if (r == UINT32_MAX) continue;
                any = true;
                const uint32_t ra = r;
                bool bonded = false;
                md_bond_iter_t it = md_bond_iter(&sys.bond, at);
                for (; md_bond_iter_has_next(&it); md_bond_iter_next(&it)) bonded |= (uint32_t)md_bond_iter_atom_index(&it) == ra;
                wrong += !bonded || ra == other;
            }
            ends += any;
        }
        both_ends += ends >= 1;     /* a C=O's O has none */
    }
    /* the sticks of the covalent class leave them out */
    {
        const uint32_t beg = L->bond_kind_offset[MD_RENDER_MOL_BOND_KIND(MD_BOND_KIND_COVALENT)], cnt = L->bond_kind_offset[MD_RENDER_MOL_BOND_KIND(MD_BOND_KIND_COVALENT) + 1] - beg;
        const uint32_t* ids = (const uint32_t*)(img + L->bond_ids);
        const uint32_t* one = (const uint32_t*)(img + L->single_bonds);
        ASSERT_TRUE(L->single_bonds != 0);
        uint32_t zeros = 0;
        for (uint32_t j = 0; j < cnt; ++j) {
            const uint32_t id = ids[beg + j], set = (one[id >> 5] >> (id & 31)) & 1u;
            wrong += (set == 0) != (is_ord[id] != 0);
            zeros += set == 0;
        }
        EXPECT_EQ(nd + nt, zeros);
    }
    printf("    1LAF: %u double, %u triple bonds; %u doubles with a plane by an end\n", nd, nt, both_ends);
    EXPECT_EQ(0, wrong);
    EXPECT_GT(both_ends, (int)nd * 9 / 10);
    /* the aromatic rings' circles: every aromatic ring of static atoms, the fused ones dashed, last */
    uint32_t n_circ = 0, n_dash = 0, n_arom = 0;
    for (size_t rg = 0; rg < md_ring_data_count(&sys.ring); ++rg) n_arom += (md_ring_flags(&sys.ring, rg) & MD_RING_FLAG_AROMATIC) != 0;
    for (uint32_t k = 0; k < MD_RENDER_MOL_RING_SIZES; ++k) {
        const uint32_t* ap = (const uint32_t*)(img + L->orders.aromatic_atoms[k]);
        const uint32_t* ai = (const uint32_t*)(img + L->orders.aromatic_id[k]);
        for (uint32_t j = 0; j < L->orders.aromatic_count[k]; ++j) {
            const md_ring_flags_t rf = md_ring_flags(&sys.ring, ai[j]);
            const bool dashed = j >= L->orders.aromatic_count[k] - L->orders.aromatic_dashed[k];
            wrong += !(rf & MD_RING_FLAG_AROMATIC) || ((rf & MD_RING_FLAG_FUSED_AROMATIC) != 0) != dashed || md_ring_size(&sys.ring, ai[j]) != 3u + k;
            const int32_t* ra = md_ring_atoms(&sys.ring, ai[j]);
            for (uint32_t l = 0; l < 3u + k; ++l) wrong += (int32_t)ap[(3u + k) * j + l] != ra[l];
        }
        n_circ += L->orders.aromatic_count[k];
        n_dash += L->orders.aromatic_dashed[k];
    }
    /* the partners: the delocalized bonds, the aromatic ones no circle covers (none here: every ring has one) */
    uint8_t* in_circle = calloc(sys.bond.count + 1, 1);
    for (uint32_t k = 0; k < MD_RENDER_MOL_RING_SIZES; ++k) {
        const uint32_t* ai = (const uint32_t*)(img + L->orders.aromatic_id[k]);
        for (uint32_t j = 0; j < L->orders.aromatic_count[k]; ++j) {
            const int32_t* ra = md_ring_atoms(&sys.ring, ai[j]);
            for (uint32_t l = 0; l < 3u + k; ++l) {
                md_bond_iter_t it = md_bond_iter(&sys.bond, (size_t)ra[l]);
                for (; md_bond_iter_has_next(&it); md_bond_iter_next(&it)) if (md_bond_iter_atom_index(&it) == ra[(l + 1) % (3u + k)]) in_circle[md_bond_iter_bond_index(&it)] = 1;
            }
        }
    }
    uint32_t n_deloc = 0, n_small = 0;
    int32_t* mol_of = molecule_of_atoms(&sys);
    for (size_t i = 0; i < sys.bond.count; ++i) {
        const uint32_t f = sys.bond.flags[i];
        if (!(md_bond_kind(f) == MD_BOND_KIND_COVALENT && (((f & MD_BOND_FLAG_DELOCALIZED) && !(f & MD_BOND_FLAG_AROMATIC)) || ((f & MD_BOND_FLAG_AROMATIC) && !in_circle[i])))) continue;
        const int32_t m = mol_of[sys.bond.pairs[i].idx[0]];      /* a small molecule's (an ion) has one too */
        const bool small = m >= 0 && sys.molecule.offset[m + 1] - sys.molecule.offset[m] <= MD_RENDER_MOL_SMALL;
        n_small += small;
        n_deloc++;
    }
    free(mol_of);
    const uint32_t* pbnd = (const uint32_t*)(img + L->orders.partner_id);
    const uint32_t* ppts = (const uint32_t*)(img + L->orders.partner_atoms);
    for (uint32_t k = 0; k < L->orders.partner_count; ++k) {
        const uint32_t f = sys.bond.flags[pbnd[k]];
        wrong += !(((f & MD_BOND_FLAG_DELOCALIZED) && !(f & MD_BOND_FLAG_AROMATIC)) || ((f & MD_BOND_FLAG_AROMATIC) && !in_circle[pbnd[k]]));
        wrong += ppts[2 * k] != (uint32_t)sys.bond.pairs[pbnd[k]].idx[0] || ppts[2 * k + 1] != (uint32_t)sys.bond.pairs[pbnd[k]].idx[1];
    }
    printf("    %u aromatic rings, %u circles (%u dashed); %u partners of %u bonds that want one (%u in small molecules)\n", n_arom, n_circ, n_dash, L->orders.partner_count, n_deloc, n_small);
    EXPECT_EQ(n_deloc, L->orders.partner_count);
    free(in_circle);
    EXPECT_EQ(0, wrong);
    EXPECT_EQ(n_arom, n_circ);
    EXPECT_GT(n_dash, 0u);
    EXPECT_LT(n_dash, n_circ);

    fix_t fx;
    if (!f_open(&fx)) { md_temp_end(temp); free(is_ord); md_system_free(&sys); UTEST_SKIP(f_skip_reason()); }
    uint64_t content[MD_RENDER_CONTENT_COUNT] = {0};
    md_render_mol_topo_t topo;
    ASSERT_TRUE(md_render_mol_write(fx.gfx, &topo, &plan, md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, plan.size).gpu, content));
    md_render_mol_ctx_t kernels = md_render_mol_create(fx.dev, NULL);
    ASSERT_TRUE(kernels != NULL);
    md_gpu_addr_t dyn = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, topo.dynamic_size).gpu;
    md_gpu_addr_t P = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, (uint64_t)sys.atom.count * 12).gpu;
    ASSERT_TRUE(md_render_mol_upload_state(fx.gfx, P, &state, content));
    md_render_mol_atoms_t atoms = {0};
    atoms.position = P;
    ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, 0, content));

    /* the sides against the reference */
    uint8_t* blk = read_block(&fx, dyn, topo.dynamic_size);
    ASSERT_TRUE(blk != NULL);
    md_render_batch_t bt[16];
    md_render_mol_rep_t rep = {0};
    rep.key = 7; rep.type = MD_RENDER_MOL_REP_LICORICE; rep.params.licorice.radius = 0.2f; rep.bond_orders = true;
    const uint32_t nb = md_render_mol_emit(&topo, dyn, &atoms, &rep, bt, 16);
    rep.bond_orders = false;
    const uint32_t nb0 = md_render_mol_emit(&topo, dyn, &atoms, &rep, NULL, 0);
    rep.bond_orders = true;
    uint32_t n_cb = 0;            /* batches of circles: per size, solid and dashed */
    for (uint32_t k = 0; k < MD_RENDER_MOL_RING_SIZES; ++k) n_cb += (topo.orders.aromatic_count[k] > topo.orders.aromatic_dashed[k]) + (topo.orders.aromatic_dashed[k] > 0);
    ASSERT_EQ(nb0 + 3u + n_cb, nb);   /* the covalent sticks, its doubles' two, the circles, the partners */
    const md_render_batch_t* bpt = &bt[3 + n_cb];
    EXPECT_EQ((uint64_t)(7 << 8 | MD_RENDER_MOL_PART_PARTNER), bpt->key);
    EXPECT_EQ(topo.orders.partner_count, bpt->count);
    EXPECT_NEAR(0.066f, bpt->style.capsules.radius, 1e-6f);
    EXPECT_NEAR(0.38f / 0.066f, bpt->style.capsules.offset, 1e-4f);
    EXPECT_NEAR(0.08f, bpt->style.capsules.trim, 1e-6f);
    EXPECT_NEAR(0.15f, bpt->style.capsules.dash_length, 1e-6f);
    EXPECT_NEAR(0.15f, bpt->style.capsules.gap_length, 1e-6f);
    const md_render_batch_t* c6 = NULL;
    for (uint32_t i = 3; i < 3 + n_cb; ++i) {
        const md_render_polygon_style_t* ps = &bt[i].style.polygons;
        EXPECT_EQ(MD_RENDER_POLYGONS, bt[i].kind);
        EXPECT_NEAR(0.066f, ps->radius, 1e-6f);
        EXPECT_NEAR(0.38f, ps->inset, 1e-6f);
        const uint32_t part = (uint32_t)(bt[i].key & 0xFF) - MD_RENDER_MOL_PART_AROMATIC;
        EXPECT_EQ(3u + part / 2u, ps->corners);
        EXPECT_EQ(part % 2u ? 2u : 0u, ps->dashes);      /* two over each bond: 1.39 bonds leave room */
        if (ps->corners == 6 && !ps->dashes) c6 = &bt[i];
    }
    ASSERT_TRUE(c6 != NULL);
    const md_render_batch_t* d0 = &bt[1];
    const md_render_batch_t* d1 = &bt[2];
    EXPECT_EQ(topo.single_bonds, bt[0].id_visible[0]);
    EXPECT_EQ((uint64_t)(7 << 8 | MD_RENDER_MOL_PART_DOUBLE), d0->key);
    EXPECT_EQ((uint64_t)(7 << 8 | (MD_RENDER_MOL_PART_DOUBLE + 1)), d1->key);
    EXPECT_NEAR(0.084f, d0->style.capsules.radius, 1e-6f);
    EXPECT_NEAR(0.58f / 0.42f, d0->style.capsules.offset, 1e-5f);
    EXPECT_NEAR(-0.58f / 0.42f, d1->style.capsules.offset, 1e-5f);
    EXPECT_EQ(d0->style.capsules.side, d1->style.capsules.side);
    EXPECT_EQ(nd, d0->count);
    const float* side = (const float*)(blk + (d0->style.capsules.side - dyn));
    double worst = 0;
    int strong = 0, first_strong = -1;
    for (uint32_t k = 0; k < nd; ++k) {
        double pa[3], pb[3], c[4][3];
        const float* xa = &state.xyz[op[2 * k]].x;
        const float* xb = &state.xyz[op[2 * k + 1]].x;
        for (int i = 0; i < 3; ++i) { pa[i] = xa[i]; pb[i] = xb[i]; }
        for (int q = 0; q < 4; ++q) {
            const uint32_t r = orf[4 * k + q];
            if (r != UINT32_MAX) for (int i = 0; i < 3; ++i) c[q][i] = (&state.xyz[r].x)[i];
        }
        double ax[3] = { pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2] }, l = sqrt(dot3(ax, ax));
        for (int i = 0; i < 3; ++i) ax[i] /= l;
        /* the kernel's basis is its own: compare the sides as lines of a weight */
        double u[3], v[3], hlp[3] = { 1, 0, 0 };
        if (fabs(ax[0]) >= 0.9) { hlp[0] = 0; hlp[1] = 1; }
        cross3(ax, hlp, u);
        double lu = sqrt(dot3(u, u));
        for (int i = 0; i < 3; ++i) u[i] /= lu;
        cross3(ax, u, v);
        double sv[2] = { 0, 0 };
        ref_end_lines(ax, pa, orf[4 * k] != UINT32_MAX ? c[0] : NULL, orf[4 * k + 1] != UINT32_MAX ? c[1] : NULL, u, v, sv);
        ref_end_lines(ax, pb, orf[4 * k + 2] != UINT32_MAX ? c[2] : NULL, orf[4 * k + 3] != UINT32_MAX ? c[3] : NULL, u, v, sv);
        const double w = sqrt(sv[0] * sv[0] + sv[1] * sv[1]), ang = 0.5 * atan2(sv[1], sv[0]);
        double ref[3];
        for (int i = 0; i < 3; ++i) ref[i] = (cos(ang) * u[i] + sin(ang) * v[i]) * w;
        const float* g = side + 4 * k;
        const double gl = sqrt((double)g[0] * g[0] + (double)g[1] * g[1] + (double)g[2] * g[2]);
        const double dt = fabs(g[0] * ref[0] + g[1] * ref[1] + g[2] * ref[2]);   /* lines: either sign */
        worst = fmax(worst, fabs(gl - w) + (w > 1e-3 ? fabs(dt / (gl * w) - 1.0) : 0.0));
        if (w > 1.0 && first_strong < 0) first_strong = (int)k;
        strong += w > 1.0;
    }
    printf("    sides: %d strong (> 1); the worst off the reference %.2e\n", strong, worst);
    EXPECT_LT(worst, 1e-4);
    EXPECT_GT(strong, (int)nd / 2);
    free(blk);

    /* The first strong one side on, the two sticks above and below: each
       picks as the bond, from its own batch */
    ASSERT_GE(first_strong, 0);
    {
        const uint32_t k = (uint32_t)first_strong;
        const float* gs = side + 0;     /* (read again below) */
        (void)gs;
        uint8_t* blk2 = read_block(&fx, dyn, topo.dynamic_size);
        const float* g = (const float*)(blk2 + (d0->style.capsules.side - dyn)) + 4 * k;
        const float* xa = &state.xyz[op[2 * k]].x;
        const float* xb = &state.xyz[op[2 * k + 1]].x;
        double x[3] = { xb[0] - xa[0], xb[1] - xa[1], xb[2] - xa[2] }, len = sqrt(dot3(x, x));
        for (int i = 0; i < 3; ++i) x[i] /= len;
        double y[3] = { g[0], g[1], g[2] }, ly = sqrt(dot3(y, y));
        for (int i = 0; i < 3; ++i) y[i] /= ly;
        double z[3];
        cross3(x, y, z);
        double eye[3];
        for (int i = 0; i < 3; ++i) eye[i] = 0.5 * (xa[i] + xb[i]) + 30.0 * z[i];
        const uint32_t W = 128, H = 128;
        cam_t c = { .ortho = true, .half_h = 1.0f, .aspect = 1.0f, .n = 0.0f, .far_ = 100.0f };
        md_render_camera_t cam = make_camera(&c);
        for (int i = 0; i < 3; ++i) { cam.view[4 * i + 0] = (float)x[i]; cam.view[4 * i + 1] = (float)y[i]; cam.view[4 * i + 2] = (float)z[i]; }
        cam.view[12] = (float)-dot3(x, eye); cam.view[13] = (float)-dot3(y, eye); cam.view[14] = (float)-dot3(z, eye);
        md_render_frame_t f = frame_of(W, H, cam, bt + 1, 2);
        memcpy(f.content, content, sizeof(content));
        ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
        const double o = 0.58 * 0.2;    /* the sticks' axes off the bond's */
        for (int m = 0; m < 2; ++m) {
            const double sy = m == 0 ? o : -o;
            const uint32_t py = (uint32_t)((0.5 - sy / c.half_h * 0.5) * H);
            md_render_pick_t pk;
            md_render_mol_hit_t hit;
            ASSERT_TRUE(f_pick(&fx, W / 2, py, &pk));
            ASSERT_TRUE(md_render_mol_decode_pick(&pk, &hit));
            EXPECT_EQ(MD_RENDER_MOL_HIT_BOND, hit.kind);
            EXPECT_EQ(ob[k], hit.bond);
            EXPECT_EQ((uint64_t)(7 << 8 | (MD_RENDER_MOL_PART_DOUBLE + m)), pk.key);
        }
        free(blk2);
    }

    /* The partners' sides and miters against md_gl's plane_side and inset_factor; then a
       carboxylate face on, the partner of one of its bonds picked along its centre line */
    {
        uint8_t* blk3 = read_block(&fx, dyn, topo.dynamic_size);
        const float*    ps   = (const float*)(blk3 + (bpt->style.capsules.side - dyn));
        const uint32_t* pref = (const uint32_t*)(img + L->orders.partner_ref);
        const uint32_t* ppt  = (const uint32_t*)(img + L->orders.partner_atoms);
        const uint32_t* pbd  = (const uint32_t*)(img + L->orders.partner_id);
        double pworst = 0;
        int carbox = -1;
        for (uint32_t k = 0; k < topo.orders.partner_count; ++k) {
            double pa[3], pb[3], ca[3] = { 0, 0, 0 }, cb[3] = { 0, 0, 0 };
            for (int i = 0; i < 3; ++i) { pa[i] = (&state.xyz[ppt[2 * k]].x)[i]; pb[i] = (&state.xyz[ppt[2 * k + 1]].x)[i]; }
            const bool ha = pref[2 * k] != UINT32_MAX, hb = pref[2 * k + 1] != UINT32_MAX;
            if (ha) for (int i = 0; i < 3; ++i) ca[i] = (&state.xyz[pref[2 * k]].x)[i] - pa[i];
            if (hb) for (int i = 0; i < 3; ++i) cb[i] = (&state.xyz[pref[2 * k + 1]].x)[i] - pb[i];
            double ax[3] = { pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2] }, l = sqrt(dot3(ax, ax));
            for (int i = 0; i < 3; ++i) ax[i] /= l;
            double nv[3] = { 0, 0, 0 }, c2 = 0, sd[3] = { 0, 0, 0 };
            if (ha) { cross3(ax, ca, nv); c2 = dot3(ca, ca); }
            if (hb) { double nb[3]; cross3(ax, cb, nb); const double sg = dot3(nv, nb) < 0 ? -1 : 1; for (int i = 0; i < 3; ++i) nv[i] += sg * nb[i]; c2 = fmax(c2, dot3(cb, cb)); }
            if (c2 > 0 && dot3(nv, nv) > 0.01 * c2) {
                cross3(nv, ax, sd);
                const double sl = sqrt(dot3(sd, sd));
                for (int i = 0; i < 3; ++i) sd[i] /= sl;
                if ((ha && dot3(sd, ca) < 0) || (!ha && hb && dot3(sd, cb) < 0)) for (int i = 0; i < 3; ++i) sd[i] = -sd[i];
            }
            double mt[2];
            for (int e = 0; e < 2; ++e) {
                const double* c = e ? cb : ca;
                if (!(e ? hb : ha)) { mt[e] = 0.57735; continue; }
                double cs = (e ? -1 : 1) * dot3(ax, c) / sqrt(dot3(c, c));
                cs = fmin(fmax(cs, -0.999), 0.999);
                mt[e] = sqrt((1 + cs) / (1 - cs));
            }
            const float* g = ps + 4 * k;
            uint32_t w;
            memcpy(&w, &g[3], 4);
            const double gm[2] = { half_to_float(w & 0xFFFFu), half_to_float(w >> 16) };
            for (int i = 0; i < 3; ++i) pworst = fmax(pworst, fabs(g[i] - sd[i]));
            for (int e = 0; e < 2; ++e) pworst = fmax(pworst, fabs(gm[e] - mt[e]) / fmax(mt[e], 1.0) * 0.1);   /* halves: 1e-3 relative */
            /* a carboxylate's (an O-C-O centre with two delocalized bonds) with a side */
            if (carbox < 0 && dot3(sd, sd) > 0.5 && ha != hb && md_atom_atomic_number(&sys.atom, ppt[2 * k + (ha ? 0 : 1)]) == 6) carbox = (int)k;
        }
        printf("    partners: the worst off the reference %.2e\n", pworst);
        EXPECT_LT(pworst, 2e-4);
        ASSERT_GE(carbox, 0);
        const uint32_t k = (uint32_t)carbox;
        const float* g = ps + 4 * k;
        double pa[3], pb[3];
        for (int i = 0; i < 3; ++i) { pa[i] = (&state.xyz[ppt[2 * k]].x)[i]; pb[i] = (&state.xyz[ppt[2 * k + 1]].x)[i]; }
        double x[3] = { pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2] }, len = sqrt(dot3(x, x));
        for (int i = 0; i < 3; ++i) x[i] /= len;
        double y[3] = { g[0], g[1], g[2] }, z[3];
        cross3(x, y, z);
        uint32_t w;
        memcpy(&w, &g[3], 4);
        const double d = 0.38, ia = d * half_to_float(w & 0xFFFFu) + 0.08 * len, ib = d * half_to_float(w >> 16) + 0.08 * len;
        double eye[3];
        for (int i = 0; i < 3; ++i) eye[i] = 0.5 * (pa[i] + pb[i]) + 30.0 * z[i];
        const uint32_t W = 160, H = 160;
        cam_t c = { .ortho = true, .half_h = 1.6f, .aspect = 1.0f, .n = 0.0f, .far_ = 100.0f };
        md_render_camera_t cam = make_camera(&c);
        for (int i = 0; i < 3; ++i) { cam.view[4 * i + 0] = (float)x[i]; cam.view[4 * i + 1] = (float)y[i]; cam.view[4 * i + 2] = (float)z[i]; }
        cam.view[12] = (float)-dot3(x, eye); cam.view[13] = (float)-dot3(y, eye); cam.view[14] = (float)-dot3(z, eye);
        md_render_frame_t f = frame_of(W, H, cam, bt, nb);
        memcpy(f.content, content, sizeof(content));
        ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
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
        /* along its centre line: from -len/2 + ia to len/2 - ib in x, d up in y */
        int on = 0, other = 0;
        for (int sI = 0; sI <= 20; ++sI) {
            const double sx = -0.5 * len + ia + (len - ia - ib) * sI / 20.0;
            const uint32_t px = (uint32_t)((sx / c.half_h * 0.5 + 0.5) * W);
            const uint32_t py = (uint32_t)((0.5 - d / c.half_h * 0.5) * H);
            md_render_pick_t pk;
            md_render_mol_hit_t hit;
            if (!f_pick(&fx, px, py, &pk) || !md_render_mol_decode_pick(&pk, &hit)) continue;
            if ((pk.key & 0xFF) == MD_RENDER_MOL_PART_PARTNER) { on += hit.kind == MD_RENDER_MOL_HIT_BOND && hit.bond == pbd[k]; other += hit.bond != pbd[k]; }
        }
        printf("    a carboxylate's partner: %d of 21 picks along it on its dashes\n", on);
        EXPECT_GT(on, 5);
        EXPECT_EQ(0, other);
        free(blk3);
    }

    /* The first solid 6-ring's circle face on: picks as the ring, the corner nearest the hit */
    {
        const uint32_t* ap = (const uint32_t*)(img + L->orders.aromatic_atoms[3]);
        const uint32_t* ai = (const uint32_t*)(img + L->orders.aromatic_id[3]);
        double q[6][3];
        for (int l = 0; l < 6; ++l) for (int i = 0; i < 3; ++i) q[l][i] = (&state.xyz[ap[l]].x)[i];
        circ_ref_t cr;
        ASSERT_TRUE(circ_ref_build(&cr, (const double (*)[3])q, 6, 0.38, 0.066, 0));
        double x[3], y[3], z[3];
        for (int i = 0; i < 3; ++i) { x[i] = cr.t[i]; y[i] = cr.b[i]; z[i] = cr.n[i]; }
        double eye[3];
        for (int i = 0; i < 3; ++i) eye[i] = cr.c[i] + 30.0 * z[i];
        const uint32_t W = 128, H = 128;
        cam_t c = { .ortho = true, .half_h = 1.5f, .aspect = 1.0f, .n = 0.0f, .far_ = 100.0f };
        md_render_camera_t cam = make_camera(&c);
        for (int i = 0; i < 3; ++i) { cam.view[4 * i + 0] = (float)x[i]; cam.view[4 * i + 1] = (float)y[i]; cam.view[4 * i + 2] = (float)z[i]; }
        cam.view[12] = (float)-dot3(x, eye); cam.view[13] = (float)-dot3(y, eye); cam.view[14] = (float)-dot3(z, eye);
        md_render_frame_t f = frame_of(W, H, cam, c6, 1);
        memcpy(f.content, content, sizeof(content));
        ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
        const double a = 0.2 * cr.ang[1];           /* nearer the first corner */
        const uint32_t px = (uint32_t)((cr.R * cos(a) / c.half_h * 0.5 + 0.5) * W);
        const uint32_t py = (uint32_t)((0.5 - cr.R * sin(a) / c.half_h * 0.5) * H);
        md_render_pick_t pk;
        md_render_mol_hit_t hit;
        ASSERT_TRUE(f_pick(&fx, px, py, &pk));
        ASSERT_TRUE(md_render_mol_decode_pick(&pk, &hit));
        EXPECT_EQ(MD_RENDER_MOL_HIT_RING, hit.kind);
        EXPECT_EQ(ai[0], hit.ring);
        EXPECT_EQ(ap[0], hit.atom);
    }

    /* A triple bond: three sticks facing the view, the middle one on the
       bond's axis, the others either side; picked as the bond */
    {
        md_allocator_i* a2 = md_vm_arena_create(MD_GIGABYTES(1));
        md_system_t ys = { .alloc = a2 };
        md_system_state_t yst = { .alloc = a2 };
        ASSERT_TRUE(md_pdb_system_init_from_str(&ys, &yst, (md_str_t){ octynal_pdb, sizeof(octynal_pdb) - 1 }, NULL));
        ASSERT_TRUE(md_system_infer(&ys, &yst, NULL));
        md_temp_scope_t t2 = md_temp_begin();
        md_render_mol_plan_t yp;
        ASSERT_TRUE(md_render_mol_plan(&yp, &ys, NULL, t2));
        EXPECT_EQ(1u, yp.layout.orders.order_count[0]);
        ASSERT_EQ(1u, yp.layout.orders.order_count[1]);
        md_render_mol_topo_t yt;
        ASSERT_TRUE(md_render_mol_write(fx.gfx, &yt, &yp, md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, yp.size).gpu, content));
        const uint32_t* yo  = (const uint32_t*)((const uint8_t*)yp.data + yp.layout.atom_order);
        const uint32_t* yop = (const uint32_t*)((const uint8_t*)yp.data + yp.layout.orders.order_atoms);
        const uint32_t* yob = (const uint32_t*)((const uint8_t*)yp.data + yp.layout.orders.order_id);
        const uint32_t  ta = yop[2], tb = yop[3], tbond = yob[1];
        (void)yo;
        EXPECT_EQ(3, md_bond_order(ys.bond.flags[tbond]));
        md_gpu_addr_t ydyn = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, yt.dynamic_size).gpu;
        md_gpu_addr_t YP = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, (uint64_t)ys.atom.count * 12).gpu;
        ASSERT_TRUE(md_render_mol_upload_state(fx.gfx, YP, &yst, content));
        md_render_mol_atoms_t yat = {0};
        yat.position = YP;
        ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &yt, ydyn, &yat, 0, content));
        md_render_mol_rep_t yr = {0};
        yr.key = 8; yr.type = MD_RENDER_MOL_REP_LICORICE; yr.params.licorice.radius = 0.2f; yr.bond_orders = true;
        md_render_batch_t yb[16];
        const uint32_t ny = md_render_mol_emit(&yt, ydyn, &yat, &yr, yb, 16);
        ASSERT_EQ(7u, ny);          /* covalent, the double's two, the triple's three, atoms: no ring */
        for (int m = 0; m < 3; ++m) {
            const md_render_batch_t* tq = &yb[3 + m];
            EXPECT_EQ((uint64_t)(8 << 8 | (MD_RENDER_MOL_PART_TRIPLE + m)), tq->key);
            EXPECT_NEAR(0.056f, tq->style.capsules.radius, 1e-6f);
            EXPECT_NEAR((m == 0 ? 1.0f : m == 1 ? 0.0f : -1.0f) * 0.72f / 0.28f, tq->style.capsules.offset, 1e-5f);
            EXPECT_EQ((md_gpu_addr_t)0, tq->style.capsules.side);
        }
        /* along x on the screen, orthographic: the sticks off by 0.144 down (+) and up (-) */
        const float* xa = &yst.xyz[ta].x;
        const float* xb = &yst.xyz[tb].x;
        double x[3] = { xb[0] - xa[0], xb[1] - xa[1], xb[2] - xa[2] }, len = sqrt(dot3(x, x));
        for (int i = 0; i < 3; ++i) x[i] /= len;
        double y[3] = { 0, 1, 0 }, z[3];
        if (fabs(x[1]) > 0.9) { y[0] = 1; y[1] = 0; }
        cross3(x, y, z);
        double lz = sqrt(dot3(z, z));
        for (int i = 0; i < 3; ++i) z[i] /= lz;
        cross3(z, x, y);
        double eye[3];
        for (int i = 0; i < 3; ++i) eye[i] = 0.5 * (xa[i] + xb[i]) + 30.0 * z[i];
        const uint32_t W = 128, H = 128;
        cam_t c = { .ortho = true, .half_h = 0.6f, .aspect = 1.0f, .n = 0.0f, .far_ = 100.0f };
        md_render_camera_t cam = make_camera(&c);
        for (int i = 0; i < 3; ++i) { cam.view[4 * i + 0] = (float)x[i]; cam.view[4 * i + 1] = (float)y[i]; cam.view[4 * i + 2] = (float)z[i]; }
        cam.view[12] = (float)-dot3(x, eye); cam.view[13] = (float)-dot3(y, eye); cam.view[14] = (float)-dot3(z, eye);
        md_render_frame_t f = frame_of(W, H, cam, yb + 3, 3);
        memcpy(f.content, content, sizeof(content));
        ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
        const double o = 0.72 * 0.2;
        for (int m = 0; m < 3; ++m) {
            const double sy = m == 0 ? -o : m == 1 ? 0.0 : o;       /* +offset: cross(axis, view) = -y */
            const uint32_t py = (uint32_t)((0.5 - sy / c.half_h * 0.5) * H);
            md_render_pick_t pk;
            md_render_mol_hit_t hit;
            ASSERT_TRUE(f_pick(&fx, W / 2, py, &pk));
            ASSERT_TRUE(md_render_mol_decode_pick(&pk, &hit));
            EXPECT_EQ(MD_RENDER_MOL_HIT_BOND, hit.kind);
            EXPECT_EQ(tbond, hit.bond);
            EXPECT_EQ((uint64_t)(8 << 8 | (MD_RENDER_MOL_PART_TRIPLE + m)), pk.key);
        }
        md_temp_end(t2);
        md_system_free(&ys);
        md_vm_arena_destroy(a2);
    }
    md_temp_end(temp);
    md_render_mol_destroy(kernels);
    free(is_ord);
    md_system_free(&sys);
    f_close(&fx);
}

/* The cartoon of coarse grained systems (viamd's glyphs): a Martini 3 RNA
   duplex, its bases the rings of ideal ones where their beads are (a ring
   of six per nucleotide, one of five more per purine), planar with bonds
   about as long as a base's; a Martini 3 membrane, its lipids' spheres (one
   per head, or a split head's parts) and tubes, a cholesterol's slab with
   its axes the update's; a head picked as a lipid, its atom the one its
   point stands for */
UTEST(render, mol_cartoon_coarse_grained) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    md_render_mol_ctx_t kernels = md_render_mol_create(fx.dev, NULL);
    ASSERT_TRUE(kernels != NULL);
    uint64_t content[MD_RENDER_CONTENT_COUNT] = {0};
    {
        md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
        md_system_t sys = { .alloc = alloc };
        md_system_state_t state = { .alloc = alloc };
        ASSERT_TRUE(load_system(&sys, &state, MD_UNITTEST_DATA_DIR "/martini/martini3_rna.gro"));
        uint32_t nuc = 0, purines = 0;
        const md_nucleic_backbone_data_t* NB = &sys.nucleic_backbone;
        for (size_t g = 0; g < NB->segment.count; ++g) {
            md_str_t nm = md_system_residue_name(&sys, (size_t)NB->segment.res_idx[g]);
            while (nm.len > 1 && (nm.ptr[nm.len - 1] == '3' || nm.ptr[nm.len - 1] == '5')) nm.len--;
            const char c = nm.ptr[nm.len - 1];
            nuc++;
            purines += c == 'A' || c == 'G';
        }
        md_temp_scope_t temp = md_temp_begin();
        md_render_mol_plan_t plan;
        ASSERT_TRUE(md_render_mol_plan(&plan, &sys, NULL, temp));
        const md_render_mol_topo_t* T = &plan.layout;
        printf("    RNA: %u nucleotides (%u purines): band %u, rings of six %u and of five %u, sticks %u; %u derived points\n",
               nuc, purines, T->nucleic_count, T->ring_count[3], T->ring_count[2], T->stick_count, T->derived_count);
        EXPECT_EQ(nuc, T->nucleic_count);
        EXPECT_EQ(nuc, T->ring_count[3]);
        EXPECT_EQ(purines, T->ring_count[2]);
        EXPECT_EQ(nuc, T->stick_count);
        EXPECT_EQ(6u * nuc + 3u * purines, T->derived_count);
        md_render_mol_topo_t topo;
        ASSERT_TRUE(md_render_mol_write(fx.gfx, &topo, &plan, md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, plan.size).gpu, content));
        uint32_t* r6 = malloc((size_t)T->ring_count[3] * 24 + 4);
        memcpy(r6, (const uint8_t*)plan.data + T->ring_points[3], (size_t)T->ring_count[3] * 24);
        const uint32_t n6 = T->ring_count[3];
        md_temp_end(temp);
        md_gpu_addr_t dyn = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, topo.dynamic_size).gpu;
        md_gpu_addr_t P = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, (uint64_t)sys.atom.count * 12).gpu;
        ASSERT_TRUE(md_render_mol_upload_state(fx.gfx, P, &state, content));
        md_render_mol_atoms_t atoms = {0};
        atoms.position = P;
        md_render_mol_frame_cell(atoms.cell, &state.unitcell);
        ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, 0, content));
        float* pts = (float*)read_block(&fx, dyn, ((size_t)topo.point_count + topo.derived_count) * 16);
        double lmin = 1e9, lmax = 0, worst_plane = 0;
        for (uint32_t r = 0; r < n6; ++r) {
            double q[6][3];
            for (int k = 0; k < 6; ++k) for (int i = 0; i < 3; ++i) q[k][i] = pts[4 * r6[6 * r + k] + i];
            circ_ref_t cr;
            circ_ref_build(&cr, (const double (*)[3])q, 6, 0.0, 0.0, 0);
            for (int k = 0; k < 6; ++k) {
                double e[3] = { q[(k + 1) % 6][0] - q[k][0], q[(k + 1) % 6][1] - q[k][1], q[(k + 1) % 6][2] - q[k][2] }, l = sqrt(dot3(e, e));
                lmin = fmin(lmin, l); lmax = fmax(lmax, l);
                double dc[3] = { q[k][0] - cr.c[0], q[k][1] - cr.c[1], q[k][2] - cr.c[2] };
                worst_plane = fmax(worst_plane, fabs(dot3(dc, cr.n)));
            }
        }
        printf("    the rings of six: bonds %.2f .. %.2f A, off their planes %.3f A at most\n", lmin, lmax, worst_plane);
        EXPECT_GT(lmin, 1.0);
        EXPECT_LT(lmax, 1.8);
        EXPECT_LT(worst_plane, 0.25);
        free(pts); free(r6);
        md_system_free(&sys);
        md_vm_arena_destroy(alloc);
    }
    {
        md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
        md_system_t sys = { .alloc = alloc };
        md_system_state_t state = { .alloc = alloc };
        ASSERT_TRUE(load_system(&sys, &state, MD_UNITTEST_DATA_DIR "/martini/martini3_membrane.gro"));
        md_temp_scope_t temp = md_temp_begin();
        md_render_mol_plan_t plan;
        ASSERT_TRUE(md_render_mol_plan(&plan, &sys, NULL, temp));
        md_render_mol_topo_t topo;
        ASSERT_TRUE(md_render_mol_write(fx.gfx, &topo, &plan, md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, plan.size).gpu, content));
        const uint32_t nsl = topo.lipid_tube_count[2];
        ASSERT_GT(nsl, 0u);
        uint32_t* slab_pts  = malloc((size_t)nsl * 4);
        uint32_t* slab_side = malloc((size_t)nsl * 4);
        memcpy(slab_pts, (const uint8_t*)plan.data + plan.layout.lipid_tubes + ((size_t)topo.lipid_tube_count[0] + topo.lipid_tube_count[1]) * 4, (size_t)nsl * 4);
        memcpy(slab_side, (const uint8_t*)plan.data + plan.layout.lipid_slab_side, (size_t)nsl * 4);
        md_temp_end(temp);
        md_gpu_addr_t dyn = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, topo.dynamic_size).gpu;
        md_gpu_addr_t P = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, (uint64_t)sys.atom.count * 12).gpu;
        ASSERT_TRUE(md_render_mol_upload_state(fx.gfx, P, &state, content));
        md_render_mol_atoms_t atoms = {0};
        atoms.position = P;
        md_render_mol_frame_cell(atoms.cell, &state.unitcell);
        ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, 0, content));
        md_render_mol_rep_t rep = {0};
        rep.key = 12; rep.type = MD_RENDER_MOL_REP_CARTOON;
        md_render_batch_t bt[16];
        const uint32_t nb = md_render_mol_emit(&topo, dyn, &atoms, &rep, bt, 16);
        ASSERT_EQ(4u, nb);          /* the spheres: always, a split head as one; the tubes: chains, slabs */
        EXPECT_EQ((uint64_t)(12 << 8 | MD_RENDER_MOL_PART_LIPID_SPHERES), bt[0].key);
        EXPECT_EQ(topo.lipid_sphere_count[1], bt[1].count);
        EXPECT_EQ(0u, (uint32_t)bt[1].points.radius_table);
        EXPECT_NEAR(1.1f, bt[2].style.splines.radius, 1e-6f);
        EXPECT_EQ((uint64_t)(12 << 8 | (MD_RENDER_MOL_PART_LIPID_TUBES + 2)), bt[3].key);
        EXPECT_EQ(topo.lipid_slab_profile, bt[3].style.splines.profile);
        rep.params.cartoon.split_heads = true;
        rep.params.cartoon.lipid_scale = 0.5f;
        md_render_batch_t bs[16];
        ASSERT_EQ(4u, md_render_mol_emit(&topo, dyn, &atoms, &rep, bs, 16));
        EXPECT_EQ(topo.lipid_sphere_count[2], bs[1].count);
        EXPECT_NEAR(0.5f, bs[1].style.spheres.radius_scale, 1e-6f);
        EXPECT_NEAR(0.55f, bs[2].style.splines.radius, 1e-6f);
        EXPECT_NEAR(0.5f, bs[3].style.splines.scale, 1e-6f);
        /* the slabs' axes: from each control point to its side */
        float* pts = (float*)read_block(&fx, dyn, ((size_t)topo.point_count + topo.derived_count) * 16);
        float* ax  = (float*)read_block(&fx, bt[3].style.splines.axis, (size_t)nsl * 16);
        double worst = 0;
        for (uint32_t j = 0; j < nsl; ++j) {
            const uint32_t p = slab_pts[j] & ~(uint32_t)(MD_RENDER_SPLINE_FIRST | MD_RENDER_SPLINE_LAST);
            for (int i = 0; i < 3; ++i) worst = fmax(worst, fabs(ax[4 * j + i] - (slab_side[j] != UINT32_MAX ? pts[4 * slab_side[j] + i] - pts[4 * p + i] : 0.0f)));
        }
        printf("    membrane: %zu lipids; spheres %u + %u (or %u parts), control points %u + %u; the slabs' axes within %.2e\n", sys.lipid.count,
               topo.lipid_sphere_count[0], topo.lipid_sphere_count[1], topo.lipid_sphere_count[2], topo.lipid_tube_count[0], nsl, worst);
        EXPECT_LT(worst, 1e-5);
        /* a head's sphere, the first, face on: a LIPID hit, its point's atom */
        uint32_t* sph = (uint32_t*)read_block(&fx, topo.lipid_spheres, 4);
        uint32_t* ord = (uint32_t*)read_block(&fx, topo.atom_order, ((size_t)topo.point_count + topo.derived_count) * 4);
        const float* c = pts + 4 * sph[0];
        const uint32_t W = 96, H = 96;
        cam_t cm = { .ortho = true, .half_h = 8.0f, .aspect = 1.0f, .n = 0.0f, .far_ = 1000.0f };
        md_render_camera_t cam = make_camera(&cm);
        cam.view[12] = -c[0]; cam.view[13] = -c[1]; cam.view[14] = -(c[2] + 200.0f);
        md_render_frame_t f = frame_of(W, H, cam, bt, 1);
        memcpy(f.content, content, sizeof(content));
        ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
        md_render_pick_t pk;
        md_render_mol_hit_t hit;
        ASSERT_TRUE(f_pick(&fx, W / 2, H / 2, &pk));
        ASSERT_TRUE(md_render_mol_decode_pick(&pk, &hit));
        EXPECT_EQ(MD_RENDER_MOL_HIT_LIPID, hit.kind);
        EXPECT_EQ(ord[sph[0]], hit.atom);
        free(pts); free(ax); free(sph); free(ord); free(slab_pts); free(slab_side);
        md_system_free(&sys);
        md_vm_arena_destroy(alloc);
    }
    md_render_mol_destroy(kernels);
    f_close(&fx);
}

/* RIBBONS (md_gl's): the backbone's splines as a flat stadium, oriented as
   the cartoon; its gaps as thick as the ribbon; drawn and picked as the
   backbone */
UTEST(render, mol_ribbons) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_mmcif_system_init_from_file(&sys, &state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/9gfu.cif"), NULL));
    ASSERT_TRUE(md_system_infer(&sys, &state, NULL));
    uint64_t content[MD_RENDER_CONTENT_COUNT] = {0};
    md_temp_scope_t temp = md_temp_begin();
    md_render_mol_plan_t plan;
    ASSERT_TRUE(md_render_mol_plan(&plan, &sys, NULL, temp));
    md_render_mol_topo_t topo;
    ASSERT_TRUE(md_render_mol_write(fx.gfx, &topo, &plan, md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, plan.size).gpu, content));
    md_temp_end(temp);
    md_render_mol_ctx_t kernels = md_render_mol_create(fx.dev, NULL);
    ASSERT_TRUE(kernels != NULL);
    md_gpu_addr_t dyn = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, topo.dynamic_size).gpu;
    md_gpu_addr_t P = md_gpu_malloc(fx.gfx, MD_GPU_MEM_DEVICE, (uint64_t)sys.atom.count * 12).gpu;
    ASSERT_TRUE(md_render_mol_upload_state(fx.gfx, P, &state, content));
    md_render_mol_atoms_t atoms = {0};
    atoms.position = P;
    ASSERT_TRUE(md_render_mol_update(kernels, fx.gfx, &topo, dyn, &atoms, 0, content));
    md_render_mol_rep_t rep = {0};
    rep.key = 13; rep.type = MD_RENDER_MOL_REP_RIBBONS; rep.chain_gaps = true;
    md_render_batch_t bt[4];
    ASSERT_EQ(2u, md_render_mol_emit(&topo, dyn, &atoms, &rep, bt, 4));
    EXPECT_EQ((uint64_t)(13 << 8 | MD_RENDER_MOL_PART_BACKBONE), bt[0].key);
    EXPECT_EQ(MD_RENDER_SPLINES, bt[0].kind);
    EXPECT_EQ(topo.backbone_count, bt[0].count);
    EXPECT_NEAR(0.1f, bt[0].style.splines.radius, 1e-6f);
    EXPECT_NEAR(0.9f, bt[0].style.splines.flat, 1e-6f);
    EXPECT_EQ((md_gpu_addr_t)0, bt[0].style.splines.profile);
    EXPECT_NEAR(0.14f, bt[1].style.capsules.radius, 1e-6f);
    rep.params.ribbons.half_width = 1.5f; rep.params.ribbons.half_thickness = 0.25f;
    ASSERT_EQ(2u, md_render_mol_emit(&topo, dyn, &atoms, &rep, bt, 4));
    EXPECT_NEAR(1.25f, bt[0].style.splines.flat, 1e-6f);
    EXPECT_NEAR(0.35f, bt[1].style.capsules.radius, 1e-6f);
    /* drawn, orthographic over all of it: a pick on a covered pixel is the backbone */
    float lo[3] = { 1e30f, 1e30f, 1e30f }, hi[3] = { -1e30f, -1e30f, -1e30f };
    for (size_t i = 0; i < sys.atom.count; ++i) for (int k = 0; k < 3; ++k) { lo[k] = fminf(lo[k], (&state.xyz[i].x)[k]); hi[k] = fmaxf(hi[k], (&state.xyz[i].x)[k]); }
    const uint32_t W = 160, H = 160;
    cam_t c = { .ortho = true, .half_h = 0.55f * fmaxf(hi[0] - lo[0], hi[1] - lo[1]), .aspect = 1.0f, .n = 0.0f, .far_ = 1000.0f };
    md_render_camera_t cam = make_camera(&c);
    cam.view[12] = -0.5f * (lo[0] + hi[0]); cam.view[13] = -0.5f * (lo[1] + hi[1]); cam.view[14] = -(hi[2] + 20.0f);
    md_render_frame_t f = frame_of(W, H, cam, bt, 1);
    memcpy(f.content, content, sizeof(content));
    ASSERT_TRUE(md_render_view_draw(fx.gfx, fx.view, &f, NULL));
    uint32_t* id = malloc(W * H * 8);
    ASSERT_TRUE(f_read(&fx, md_render_view_targets(fx.view).id, id, W * H * 8));
    uint32_t covered = 0, first = UINT32_MAX;
    for (uint32_t i = 0; i < W * H; ++i) if (id[2 * i + 1]) { covered++; if (first == UINT32_MAX) first = i; }
    printf("    9gfu ribbons: %u pixels of %u covered\n", covered, W * H);
    EXPECT_GT(covered, W * H / 20);
    ASSERT_NE(UINT32_MAX, first);
    md_render_pick_t pk;
    md_render_mol_hit_t hit;
    ASSERT_TRUE(f_pick(&fx, first % W, first / W, &pk));
    ASSERT_TRUE(md_render_mol_decode_pick(&pk, &hit));
    EXPECT_EQ(MD_RENDER_MOL_HIT_BACKBONE, hit.kind);
    EXPECT_LT(hit.segment, (uint32_t)sys.protein_backbone.segment.count);
    free(id);
    md_render_mol_destroy(kernels);
    md_system_free(&sys);
    f_close(&fx);
}

/* ---- One submission against each pass synced -----------------------------------
   A frame recorded in one go must draw what the same frame draws with the
   stream synced after each of its passes. On the HD 530 (Vulkan 1.3.215) it
   did not: synced pass by pass the rasters drew, in one go they left the
   targets empty, until md_gpu's implicit barrier named its stages and
   accesses (md_vk_stream_full_masks). Each variant draws on a view of its own (its first frame,
   as the tests' frames are); the lines it prints say where an ordering is
   missing: a flush (a submission boundary, no wait) or a sync after one pass
   alone that repairs the frame names the boundary. */

typedef struct order_probe_t {
    md_gpu_stream_t s;
    int             mode;           /* 0 nothing, 1 sync, 2 flush (submit, no wait) */
    const char*     at;             /* NULL: after every pass, else after this one */
    const char*     seen[32];
    int             n;
} order_probe_t;

static void order_hook(void* user, const char* pass) {
    order_probe_t* p = user;
    if (p->n < 32) p->seen[p->n++] = pass;
    if (p->at && strcmp(p->at, pass) != 0) return;
    if (p->mode == 1) md_gpu_stream_sync(p->s);
    else if (p->mode == 2) md_gpu_stream_flush(p->s);
}

typedef struct order_result_t { int ids, depth, color; uint32_t hw; bool ok; } order_result_t;

static order_result_t order_draw(fix_t* fx, const md_render_frame_t* f, order_probe_t* p) {
    order_result_t o = {0};
    const uint32_t W = f->width, H = f->height;
    md_render_view_t v = md_render_view_create(fx->r, "ordering view");
    if (!v) return o;
    p->s = fx->gfx;
    p->n = 0;
    md_render_debug_pass_hook(order_hook, p);
    md_render_status_t st;
    o.ok = md_render_view_draw(fx->gfx, v, f, &st);
    md_render_debug_pass_hook(NULL, NULL);
    md_render_targets_t t = md_render_view_targets(v);
    uint32_t* id  = malloc(W * H * 8);
    float*    dp  = malloc(W * H * 4);
    uint32_t* col = malloc(W * H * 4);
    if (o.ok && f_read(fx, t.id, id, W * H * 8) && f_read(fx, t.depth, dp, W * H * 4) && f_read(fx, t.color, col, W * H * 4)) {
        for (uint32_t i = 0; i < W * H; ++i) {
            o.ids   += (id[2 * i + 1] >> 20) != 0;
            o.depth += dp[i] > 0.0f;
            o.color += col[i] != 0;
        }
    } else {
        o.ok = false;
    }
    md_render_stats_t rs;
    if (md_render_view_stats(v, &rs)) o.hw = rs.hardware[0] + rs.hardware[1];
    free(id); free(dp); free(col);
    md_render_view_destroy(fx->gfx, v);
    md_gpu_stream_sync(fx->gfx);
    return o;
}

static void order_print(const char* what, order_result_t r) {
    printf("    %-34s ids %5d px, depth %5d px, colour %5d px, hardware %u%s\n", what, r.ids, r.depth, r.color, r.hw, r.ok ? "" : " (draw or readback failed)");
}

UTEST(render, one_submission_draws_what_synced_passes_draw) {
    fix_t fx;
    if (!f_open(&fx)) UTEST_SKIP(f_skip_reason());
    const uint32_t W = 96, H = 64;
    cam_t c = { .tan_half = 0.4f, .aspect = (float)W / H, .n = 0.1f };
    md_render_camera_t cam = make_camera(&c);
    const float pos[9]    = { 0, 0, -10,  2.2f, 0.8f, -12,  -1.5f, -0.9f, -8 };
    const float radius[3] = { 1.5f, 2.0f, 0.7f };
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
    md_gpu_stream_sync(fx.gfx);
    md_render_frame_t f = frame_of(W, H, cam, &b, 1);

    order_probe_t p = {0};
    p.mode = 1;
    const order_result_t synced = order_draw(&fx, &f, &p);
    order_print("each pass synced:", synced);
    const char* passes[32];
    const int   npass = p.n;
    for (int i = 0; i < npass; ++i) passes[i] = p.seen[i];

    p.mode = 0;
    const order_result_t one = order_draw(&fx, &f, &p);
    order_print("one submission:", one);

    if (one.ids != synced.ids || one.depth != synced.depth || one.color != synced.color) {
        p.mode = 2;
        order_print("each pass flushed (no wait):", order_draw(&fx, &f, &p));
        for (int i = 0; i < npass; ++i) {
            char what[64];
            p.at = passes[i];
            p.mode = 1;
            snprintf(what, sizeof(what), "synced after '%s' alone:", passes[i]);
            order_print(what, order_draw(&fx, &f, &p));
            p.mode = 2;
            snprintf(what, sizeof(what), "flushed after '%s' alone:", passes[i]);
            order_print(what, order_draw(&fx, &f, &p));
        }
        p.at = NULL;
    }
    EXPECT_TRUE(synced.ok && one.ok);
    EXPECT_GT(synced.ids, (int)(W * H / 10));
    EXPECT_EQ(synced.ids, one.ids);
    EXPECT_EQ(synced.depth, one.depth);
    EXPECT_EQ(synced.color, one.color);
    md_gpu_free(fx.gfx, b.points.position);
    md_gpu_free(fx.gfx, b.points.radius);
    f_close(&fx);
}

#endif
