#include "md_post.h"

#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_log.h>

#include <math.h>

#include "md_post_shaders.inl"

/* md_post.slang's Args */
typedef struct gpu_post_args_t {
    md_gpu_float4        bg, tm, pa, pb, pc, jit;
    md_gpu_float4        vis0, vis1, hid0, hid1, out0, out1;
    md_gpu_float4        g0, g1, g2, g3, g4, g5, g6;
    md_gpu_float4        mk;
    uint32_t             in_w, in_h, out_w, out_h;
    uint32_t             has_prev, img_w, img_h, _p0;
    md_gpu_sampled_tex_t src, prev, depth, motion, mask;
    md_gpu_sampler_t     bilinear;
    md_gpu_storage_tex_t dst;
} gpu_post_args_t;
_Static_assert(sizeof(gpu_post_args_t) == 416, "Args mirror");

struct md_post {
    md_gpu_device_t  dev;
    md_gpu_kernel_t  k_display, k_resolve, k_present;
    md_gpu_sampler_t bilinear;
};

md_post_ctx_t md_post_create(md_gpu_device_t device) {
    if (!device) {
        MD_LOG_ERROR("md_post: no device");
        return NULL;
    }
    md_post_ctx_t p = md_alloc(md_get_heap_allocator(), sizeof(*p));
    if (!p) return NULL;
    md_memset(p, 0, sizeof(*p));
    p->dev = device;
    md_gpu_kernel_desc_t kd;
    kd = md_shader_md_post_display_kernel(); p->k_display = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_post_resolve_kernel(); p->k_resolve = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_post_present_kernel(); p->k_present = md_gpu_kernel_create(device, &kd);
    md_gpu_sampler_desc_t sd = {0};
    sd.min_filter = sd.mag_filter = MD_GPU_FILTER_LINEAR;
    p->bilinear = md_gpu_sampler(device, &sd);
    if (!p->k_display || !p->k_resolve || !p->k_present || !p->bilinear.handle) {
        MD_LOG_ERROR("md_post: kernel creation failed: %s", md_gpu_last_error());
        md_post_destroy(p);
        return NULL;
    }
    return p;
}

void md_post_destroy(md_post_ctx_t p) {
    if (!p) return;
    if (p->k_display) md_gpu_kernel_destroy(p->k_display);
    if (p->k_resolve) md_gpu_kernel_destroy(p->k_resolve);
    if (p->k_present) md_gpu_kernel_destroy(p->k_present);
    md_free(md_get_heap_allocator(), p, sizeof(*p));
}

static md_gpu_float4 f4(float x, float y, float z, float w) {
    md_gpu_float4 v = { x, y, z, w };
    return v;
}

static md_gpu_float4 f4v(const float v[4]) { return f4(v[0], v[1], v[2], v[3]); }

/* A 2D texture of the format and usage, and (w, h unless 0) of the size */
static const md_gpu_texture_desc_t* tex_ok(md_gpu_texture_t t, md_gpu_format_t format, md_gpu_tex_usage_t usage, uint32_t w, uint32_t h) {
    const md_gpu_texture_desc_t* d = t ? md_gpu_texture_desc(t) : NULL;
    if (!d || d->type != MD_GPU_TEX_2D || d->format != format || (d->usage & usage) != usage) return NULL;
    if (w && (d->width != w || d->height != h)) return NULL;
    return d;
}

static bool launch(md_gpu_stream_t s, md_gpu_kernel_t k, uint32_t w, uint32_t h, const gpu_post_args_t* a) {
    const md_gpu_ordering_t prev = md_gpu_stream_set_ordering(s, MD_GPU_ORDER_IMPLICIT);
    const bool ok = md_gpu_launch(s, k, md_gpu_grid_for(k, w, h, 1), a, sizeof(*a));
    md_gpu_stream_set_ordering(s, prev);
    if (!ok) MD_LOG_ERROR("md_post: %s", md_gpu_last_error());
    return ok;
}

static float srgb_decode(float x) {
    x = x < 0.0f ? 0.0f : x > 1.0f ? 1.0f : x;
    return x <= 0.04045f ? x / 12.92f : powf((x + 0.055f) / 1.055f, 2.4f);
}

bool md_post_display(md_gpu_stream_t s, md_post_ctx_t p, const md_post_display_t* d) {
    if (!s || !p || !d) return false;
    const md_gpu_texture_desc_t* sd = tex_ok(d->scene, MD_GPU_FORMAT_RGBA16_FLOAT, MD_GPU_TEX_SAMPLED, 0, 0);
    if (!sd) {
        MD_LOG_ERROR("md_post_display: the scene is to be a SAMPLED 2D texture of RGBA16_FLOAT (md_compose's target)");
        return false;
    }
    if (!tex_ok(d->target, MD_GPU_FORMAT_RGBA16_FLOAT, MD_GPU_TEX_STORAGE, sd->width, sd->height)) {
        MD_LOG_ERROR("md_post_display: the target is to be a STORAGE 2D texture of RGBA16_FLOAT, %u x %u", sd->width, sd->height);
        return false;
    }
    const bool  curved = d->curve != MD_POST_CURVE_NONE;
    const float exposure = curved ? (d->exposure > 1.0e-6f ? d->exposure : 1.0e-6f) : 1.0f;
    const float gamma = curved ? (d->gamma < 0.1f ? 0.1f : d->gamma > 10.0f ? 10.0f : d->gamma) : 2.2f;
    /* the background's display referred linear value: what the encoding (with
       its gamma) shows as the colour picked */
    const float k = gamma / 2.2f;
    float bg[3];
    for (int i = 0; i < 3; ++i) {
        const float c = d->background[i] < 0.0f ? 0.0f : d->background[i] > 1.0f ? 1.0f : d->background[i];
        bg[i] = srgb_decode(powf(c, k));
    }
    gpu_post_args_t a;
    md_memset(&a, 0, sizeof(a));
    a.bg   = f4(bg[0], bg[1], bg[2], exposure);
    a.tm   = f4((float)d->curve, gamma, 0.0f, 0.0f);
    a.in_w = sd->width;
    a.in_h = sd->height;
    a.src  = md_gpu_texture_sampled(d->scene);
    a.dst  = md_gpu_texture_storage(d->target, 0);
    return launch(s, p->k_display, a.in_w, a.in_h, &a);
}

bool md_post_resolve(md_gpu_stream_t s, md_post_ctx_t p, const md_post_resolve_t* d) {
    if (!s || !p || !d || !d->view) return false;
    const md_render_targets_t t = md_render_view_targets(d->view);
    const md_gpu_texture_desc_t* vd = t.depth ? md_gpu_texture_desc(t.depth) : NULL;
    if (!vd) {
        MD_LOG_ERROR("md_post_resolve: the view has drawn nothing yet");
        return false;
    }
    if (!tex_ok(d->frame, MD_GPU_FORMAT_RGBA16_FLOAT, MD_GPU_TEX_SAMPLED, vd->width, vd->height)) {
        MD_LOG_ERROR("md_post_resolve: the frame is to be a SAMPLED 2D texture of RGBA16_FLOAT of the view's size, %u x %u", vd->width, vd->height);
        return false;
    }
    const md_gpu_texture_desc_t* hd = tex_ok(d->history, MD_GPU_FORMAT_RGBA16_FLOAT, MD_GPU_TEX_STORAGE, 0, 0);
    if (!hd || hd->width < vd->width || hd->height < vd->height) {
        MD_LOG_ERROR("md_post_resolve: the history is to be a STORAGE 2D texture of RGBA16_FLOAT, at least the view's size");
        return false;
    }
    if (d->history_prev && !tex_ok(d->history_prev, MD_GPU_FORMAT_RGBA16_FLOAT, MD_GPU_TEX_SAMPLED, hd->width, hd->height)) {
        MD_LOG_ERROR("md_post_resolve: the previous history is to be a SAMPLED 2D texture of RGBA16_FLOAT, of the history's size");
        return false;
    }
    if (d->history_prev && d->history_prev == d->history) {
        MD_LOG_ERROR("md_post_resolve: the history is to be another texture than the previous one");
        return false;
    }
    const float* P = d->camera.proj;
    const bool ortho = P[3] == 0.0f && P[7] == 0.0f && P[11] == 0.0f && P[15] == 1.0f;
    gpu_post_args_t a;
    md_memset(&a, 0, sizeof(a));
    a.tm       = f4(0.0f, 2.2f, d->feedback_min, d->feedback_max);
    a.pa       = f4(P[0], P[5], P[8], P[9]);
    a.pb       = f4(P[12], P[13], P[10], P[14]);
    a.pc       = f4(P[11], P[15], 2.0f / (P[5] * (float)vd->height), ortho ? 1.0f : 0.0f);
    a.jit      = f4(d->jitter[0], d->jitter[1], 0.0f, 0.0f);
    a.in_w     = vd->width;
    a.in_h     = vd->height;
    a.out_w    = hd->width;
    a.out_h    = hd->height;
    a.has_prev = d->history_prev ? 1u : 0u;
    a.src      = md_gpu_texture_sampled(d->frame);
    a.prev     = md_gpu_texture_sampled(d->history_prev ? d->history_prev : d->frame);
    a.depth    = md_gpu_texture_sampled(t.depth);
    a.motion   = md_gpu_texture_sampled(t.motion);
    a.bilinear = p->bilinear;
    a.dst      = md_gpu_texture_storage(d->history, 0);
    return launch(s, p->k_resolve, a.out_w, a.out_h, &a);
}

bool md_post_present(md_gpu_stream_t s, md_post_ctx_t p, const md_post_present_t* d) {
    if (!s || !p || !d) return false;
    const md_gpu_texture_desc_t* id = tex_ok(d->image, MD_GPU_FORMAT_RGBA16_FLOAT, MD_GPU_TEX_SAMPLED, 0, 0);
    if (!id) {
        MD_LOG_ERROR("md_post_present: the image is to be a SAMPLED 2D texture of RGBA16_FLOAT");
        return false;
    }
    const md_gpu_texture_desc_t* td = tex_ok(d->target, MD_GPU_FORMAT_RGBA8_UNORM, MD_GPU_TEX_STORAGE, 0, 0);
    if (!td || td->width < id->width || td->height < id->height) {
        MD_LOG_ERROR("md_post_present: the target is to be a STORAGE 2D texture of RGBA8_UNORM, at least the image's size");
        return false;
    }
    gpu_post_args_t a;
    md_memset(&a, 0, sizeof(a));
    a.out_w    = td->width;
    a.out_h    = td->height;
    a.img_w    = id->width;
    a.img_h    = id->height;
    a.src      = md_gpu_texture_sampled(d->image);
    a.bilinear = p->bilinear;
    a.dst      = md_gpu_texture_storage(d->target, 0);
    a.mk       = f4(0.0f, 1.0f, (float)(d->frame % 1000u) * 0.01f, 0.0f);
    const md_render_targets_t t = d->marks_view && d->marks ? md_render_view_targets(d->marks_view) : (md_render_targets_t){0};
    const md_gpu_texture_desc_t* md = t.mask ? md_gpu_texture_desc(t.mask) : NULL;
    a.mask = md_gpu_texture_sampled(md ? t.mask : d->image);
    if (md) {
        const md_post_marks_t* m = d->marks;
        a.in_w = md->width;
        a.in_h = md->height;
        a.vis0 = f4v(m->visible[0]); a.vis1 = f4v(m->visible[1]);
        a.hid0 = f4v(m->hidden[0]);  a.hid1 = f4v(m->hidden[1]);
        a.out0 = f4v(m->outline[0]); a.out1 = f4v(m->outline[1]);
        md_gpu_float4* g[7] = { &a.g0, &a.g1, &a.g2, &a.g3, &a.g4, &a.g5, &a.g6 };
        for (int i = 0; i < 7; ++i) *g[i] = f4(m->groups[i][0], m->groups[i][1], m->groups[i][2], 0.0f);
        const float h = m->highlight < 0.0f ? 0.0f : m->highlight > 1.0f ? 1.0f : m->highlight;
        a.mk = f4(h, (float)(m->outline_px ? m->outline_px : 1u), a.mk.z, 1.0f);
    }
    return launch(s, p->k_present, a.out_w, a.out_h, &a);
}
