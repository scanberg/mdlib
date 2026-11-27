#include "md_post.h"

#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_log.h>

#include <math.h>

#include "md_post_shaders.inl"
#include "md_post_overlay_shaders.inl"
#include "md_post_ui_shaders.inl"

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
    uint32_t             _p1, _p2;
    md_gpu_float4        mb;            /* the motion blur: scale, the noise's seed, the tiles across, down */
    md_gpu_addr_t        tmax, tiles;
    md_gpu_storage_tex_t dst2;
} gpu_post_args_t;
_Static_assert(sizeof(gpu_post_args_t) == 464, "Args mirror");

enum { VEL_TILE = 8 };  /* md_post.slang's */

/* md_post_overlay.slang's Args and Item */
typedef struct gpu_overlay_args_t {
    md_gpu_float4        v0, v1, v2, v3, q0, q1, q2, q3, sz, nf;
    md_gpu_addr_t        items;
    uint32_t             count, has_depth;
    md_gpu_sampled_tex_t depth;
    uint32_t             _p0, _p1;
} gpu_overlay_args_t;
_Static_assert(sizeof(gpu_overlay_args_t) == 192, "Overlay args mirror");
_Static_assert(sizeof(md_post_overlay_item_t) == 64, "Overlay item mirror");

/* md_post_ui.slang's Args */
typedef struct gpu_ui_args_t {
    md_gpu_float4        xf;
    md_gpu_addr_t        vertices;
    md_gpu_sampled_tex_t tex;
    md_gpu_sampler_t     smp;
    uint64_t             _p0;
} gpu_ui_args_t;
_Static_assert(sizeof(gpu_ui_args_t) == 48, "UI args mirror");
_Static_assert(sizeof(md_post_ui_vertex_t) == 20, "UI vertex: ImDrawVert's layout");

/* The UI's pipelines, by the target's format */
static const md_gpu_format_t ui_formats[2] = { MD_GPU_FORMAT_RGBA8_UNORM, MD_GPU_FORMAT_BGRA8_UNORM };

struct md_post {
    md_gpu_device_t   dev;
    md_gpu_kernel_t   k_display, k_resolve, k_present, k_tilemax, k_neighbormax;
    md_gpu_pipeline_t overlay;
    md_gpu_pipeline_t ui[2][2];    /* ui_formats', straight and premultiplied */
    md_gpu_sampler_t  bilinear;
    md_gpu_texture_t  white;       /* 1 x 1, the UI's texture where none is given: made at its first call */
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
    kd = md_shader_md_post_vel_tilemax_kernel();     p->k_tilemax     = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_post_vel_neighbormax_kernel(); p->k_neighbormax = md_gpu_kernel_create(device, &kd);
    md_gpu_sampler_desc_t sd = {0};
    sd.min_filter = sd.mag_filter = MD_GPU_FILTER_LINEAR;
    p->bilinear = md_gpu_sampler(device, &sd);
    /* The overlay's: premultiplied over the target */
    md_gpu_pipeline_desc_t pd = {0};
    pd.vertex   = md_shader_md_post_overlay_overlay_vs_shader();
    pd.fragment = md_shader_md_post_overlay_overlay_fs_shader();
    pd.topology = MD_GPU_TOPOLOGY_TRIANGLES;
    pd.color_count = 1;
    pd.color[0].format          = MD_GPU_FORMAT_RGBA8_UNORM;
    pd.color[0].blend.enable    = true;
    pd.color[0].blend.src_color = MD_GPU_BLEND_ONE;
    pd.color[0].blend.dst_color = MD_GPU_BLEND_ONE_MINUS_SRC_ALPHA;
    pd.color[0].blend.src_alpha = MD_GPU_BLEND_ONE;
    pd.color[0].blend.dst_alpha = MD_GPU_BLEND_ONE_MINUS_SRC_ALPHA;
    pd.label = "md_post.overlay";
    p->overlay = md_gpu_pipeline_create(device, &pd);
    /* The UI's: straight alpha over the target */
    md_memset(&pd, 0, sizeof(pd));
    pd.vertex   = md_shader_md_post_ui_ui_vs_shader();
    pd.fragment = md_shader_md_post_ui_ui_fs_shader();
    pd.topology = MD_GPU_TOPOLOGY_TRIANGLES;
    pd.color_count = 1;
    pd.color[0].blend.enable    = true;
    pd.color[0].blend.dst_color = MD_GPU_BLEND_ONE_MINUS_SRC_ALPHA;
    pd.color[0].blend.src_alpha = MD_GPU_BLEND_ONE;
    pd.color[0].blend.dst_alpha = MD_GPU_BLEND_ONE_MINUS_SRC_ALPHA;
    pd.label = "md_post.ui";
    bool ui_ok = true;
    for (int i = 0; i < 2; ++i) {
        for (int m = 0; m < 2; ++m) {
            pd.color[0].format = ui_formats[i];
            pd.color[0].blend.src_color = m ? MD_GPU_BLEND_ONE : MD_GPU_BLEND_SRC_ALPHA;
            p->ui[i][m] = md_gpu_pipeline_create(device, &pd);
            ui_ok = ui_ok && p->ui[i][m];
        }
    }
    if (!p->k_display || !p->k_resolve || !p->k_present || !p->k_tilemax || !p->k_neighbormax || !p->overlay || !ui_ok || !p->bilinear.handle) {
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
    if (p->k_tilemax) md_gpu_kernel_destroy(p->k_tilemax);
    if (p->k_neighbormax) md_gpu_kernel_destroy(p->k_neighbormax);
    if (p->overlay)   md_gpu_pipeline_destroy(p->overlay);
    for (int i = 0; i < 2; ++i) {
        for (int m = 0; m < 2; ++m) {
            if (p->ui[i][m]) md_gpu_pipeline_destroy(p->ui[i][m]);
        }
    }
    if (p->white)     md_gpu_texture_destroy(p->white);
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
    const bool blur = d->motion_blur > 0.0f;
    if (blur && (!tex_ok(d->shown, MD_GPU_FORMAT_RGBA16_FLOAT, MD_GPU_TEX_STORAGE, hd->width, hd->height) || d->shown == d->history)) {
        MD_LOG_ERROR("md_post_resolve: with the motion blur, the image shown is to be another STORAGE 2D texture of RGBA16_FLOAT, of the history's size");
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
    a.dst2     = md_gpu_texture_storage(blur ? d->shown : d->history, 0);
    if (!blur) return launch(s, p->k_resolve, a.out_w, a.out_h, &a);

    /* the motion reduced to tiles, then to the fastest about each, for the
       resolve to blur along */
    const uint32_t tw = (vd->width + VEL_TILE - 1) / VEL_TILE, th = (vd->height + VEL_TILE - 1) / VEL_TILE;
    a.mb = f4(d->motion_blur, (float)(d->seed % 1000u) * 0.01f, (float)tw, (float)th);
    md_gpu_temp_t scope = md_gpu_temp_begin(s);
    a.tmax  = md_gpu_temp_alloc(scope, MD_GPU_MEM_DEVICE, (size_t)tw * th * 8u).gpu;
    a.tiles = md_gpu_temp_alloc(scope, MD_GPU_MEM_DEVICE, (size_t)tw * th * 8u).gpu;
    bool ok = a.tmax && a.tiles && launch(s, p->k_tilemax, tw, th, &a) && launch(s, p->k_neighbormax, tw, th, &a) &&
              launch(s, p->k_resolve, a.out_w, a.out_h, &a);
    md_gpu_temp_end(scope);
    return ok;
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

bool md_post_overlay(md_gpu_stream_t s, md_post_ctx_t p, const md_post_overlay_t* d) {
    if (!s || !p || !d) return false;
    const md_gpu_texture_desc_t* td = tex_ok(d->target, MD_GPU_FORMAT_RGBA8_UNORM, MD_GPU_TEX_RENDER_TARGET, 0, 0);
    if (!td) {
        MD_LOG_ERROR("md_post_overlay: the target is to be a RENDER_TARGET 2D texture of RGBA8_UNORM");
        return false;
    }
    const md_gpu_texture_desc_t* dd = d->depth ? tex_ok(d->depth, MD_GPU_FORMAT_D32_FLOAT, MD_GPU_TEX_SAMPLED, 0, 0) : NULL;
    if (d->depth && !dd) {
        MD_LOG_ERROR("md_post_overlay: the depth is to be a SAMPLED 2D texture of D32_FLOAT");
        return false;
    }
    if (d->count == 0) return true;
    if (!d->items) {
        MD_LOG_ERROR("md_post_overlay: %u items and none given", d->count);
        return false;
    }
    gpu_overlay_args_t a;
    md_memset(&a, 0, sizeof(a));
    const float* V = d->camera.view;
    const float* P = d->camera.proj;
    a.v0 = f4v(V + 0);  a.v1 = f4v(V + 4);  a.v2 = f4v(V + 8);  a.v3 = f4v(V + 12);
    a.q0 = f4v(P + 0);  a.q1 = f4v(P + 4);  a.q2 = f4v(P + 8);  a.q3 = f4v(P + 12);
    a.sz = f4((float)td->width, (float)td->height, dd ? (float)dd->width : 1.0f, dd ? (float)dd->height : 1.0f);
    /* The near plane, where the depth is 1: d (P23 z + P33) = P22 z + P32 */
    const bool ortho = P[11] == 0.0f && P[15] == 1.0f;
    const float den = P[11] - P[10];
    a.nf = f4(!ortho && den != 0.0f ? -(P[14] - P[15]) / den : 0.0f, ortho ? 1.0f : 0.0f, 0.0f, 0.0f);
    a.count = d->count;
    a.has_depth = dd ? 1u : 0u;
    if (dd) a.depth = md_gpu_texture_sampled(d->depth);    /* else none: never read */

    md_gpu_temp_t scope = md_gpu_temp_begin(s);
    a.items = md_gpu_temp_upload(scope, d->items, (size_t)d->count * sizeof(md_post_overlay_item_t));
    bool ok = a.items != 0;
    if (ok) {
        md_gpu_render_desc_t rd = {0};
        rd.color_count = 1;
        rd.color[0].texture = d->target;
        rd.color[0].load  = MD_GPU_LOAD;
        rd.color[0].store = MD_GPU_STORE;
        rd.label = "md_post.overlay";
        ok = md_gpu_render_begin(s, &rd);
        if (ok) {
            ok = md_gpu_draw(s, p->overlay, 6u * d->count, 1, &a, sizeof(a));
            ok = md_gpu_render_end(s) && ok;
        }
    }
    md_gpu_temp_end(scope);
    if (!ok) MD_LOG_ERROR("md_post_overlay: %s", md_gpu_last_error());
    return ok;
}

bool md_post_ui(md_gpu_stream_t s, md_post_ctx_t p, const md_post_ui_t* d) {
    if (!s || !p || !d) return false;
    const md_gpu_texture_desc_t* td = d->target ? md_gpu_texture_desc(d->target) : NULL;
    int f = -1;
    for (int i = 0; td && i < 2; ++i) {
        if (td->format == ui_formats[i]) f = i;
    }
    if (!td || td->type != MD_GPU_TEX_2D || !(td->usage & MD_GPU_TEX_RENDER_TARGET) || f < 0) {
        MD_LOG_ERROR("md_post_ui: the target is to be a RENDER_TARGET 2D texture of RGBA8_UNORM or BGRA8_UNORM");
        return false;
    }
    if (d->list_count && !d->lists) {
        MD_LOG_ERROR("md_post_ui: %u lists and none given", d->list_count);
        return false;
    }
    const size_t isz = d->index_type == MD_GPU_INDEX_U16 ? 2u : 4u;
    for (uint32_t l = 0; l < d->list_count; ++l) {
        const md_post_ui_list_t* L = &d->lists[l];
        if ((L->vertex_count && !L->vertices) || (L->index_count && !L->indices) || (L->cmd_count && !L->cmds)) {
            MD_LOG_ERROR("md_post_ui: list %u counts what it does not give", l);
            return false;
        }
        for (uint32_t c = 0; c < L->cmd_count; ++c) {
            const md_post_ui_cmd_t* C = &L->cmds[c];
            if ((uint64_t)C->first_index + C->index_count > L->index_count) {
                MD_LOG_ERROR("md_post_ui: command %u of list %u reads indices %u .. %u of %u", c, l, C->first_index,
                             C->first_index + C->index_count, L->index_count);
                return false;
            }
        }
    }
    if (!p->white) {
        md_gpu_texture_desc_t wd = {0};
        wd.type   = MD_GPU_TEX_2D;
        wd.format = MD_GPU_FORMAT_RGBA8_UNORM;
        wd.usage  = MD_GPU_TEX_SAMPLED;
        wd.width  = wd.height = 1;
        wd.label  = "md_post.ui.white";
        const uint32_t one = 0xFFFFFFFFu;
        p->white = md_gpu_texture_create(s, &wd);
        if (!p->white || !md_gpu_upload_texture(s, p->white, NULL, &one, sizeof(one))) {
            MD_LOG_ERROR("md_post_ui: %s", md_gpu_last_error());
            return false;
        }
    }

    const float sx = d->scale[0] > 0.0f ? d->scale[0] : 1.0f, sy = d->scale[1] > 0.0f ? d->scale[1] : 1.0f;
    const float W = (float)td->width, H = (float)td->height;
    gpu_ui_args_t a;
    md_memset(&a, 0, sizeof(a));
    /* clip x = ((x - origin) scale) / W 2 - 1, y = 1 - ((y - origin) scale) / H 2 */
    a.xf = f4(2.0f * sx / W, -2.0f * sy / H, -2.0f * d->origin[0] * sx / W - 1.0f, 2.0f * d->origin[1] * sy / H + 1.0f);
    const md_gpu_sampled_tex_t white = md_gpu_texture_sampled(p->white);

    md_gpu_temp_t scope = md_gpu_temp_begin(s);
    /* Each list's vertices and indices, uploaded before the pass */
    md_gpu_addr_t* va = md_alloc(md_get_heap_allocator(), sizeof(md_gpu_addr_t) * 2 * (d->list_count ? d->list_count : 1));
    md_gpu_addr_t* ia = va + (d->list_count ? d->list_count : 1);
    bool ok = va != NULL;
    for (uint32_t l = 0; ok && l < d->list_count; ++l) {
        const md_post_ui_list_t* L = &d->lists[l];
        va[l] = L->vertex_count ? md_gpu_temp_upload(scope, L->vertices, (size_t)L->vertex_count * sizeof(md_post_ui_vertex_t)) : 0;
        ia[l] = L->index_count  ? md_gpu_temp_upload(scope, L->indices, (size_t)L->index_count * isz) : 0;
        ok = (va[l] || !L->vertex_count) && (ia[l] || !L->index_count);
    }
    if (ok) {
        md_gpu_render_desc_t rd = {0};
        rd.color_count = 1;
        rd.color[0].texture = d->target;
        rd.color[0].load  = d->clear ? MD_GPU_LOAD_CLEAR : MD_GPU_LOAD;
        rd.color[0].store = MD_GPU_STORE;
        for (int i = 0; i < 4; ++i) rd.color[0].clear.f32[i] = d->clear_color[i];
        rd.label = "md_post.ui";
        ok = md_gpu_render_begin(s, &rd);
        for (uint32_t l = 0; ok && l < d->list_count; ++l) {
            const md_post_ui_list_t* L = &d->lists[l];
            for (uint32_t c = 0; ok && c < L->cmd_count; ++c) {
                const md_post_ui_cmd_t* C = &L->cmds[c];
                if (!C->index_count) continue;
                /* The clip rectangle in the target's pixels */
                float x0 = (C->clip[0] - d->origin[0]) * sx, y0 = (C->clip[1] - d->origin[1]) * sy;
                float x1 = (C->clip[2] - d->origin[0]) * sx, y1 = (C->clip[3] - d->origin[1]) * sy;
                x0 = x0 < 0.0f ? 0.0f : x0;
                y0 = y0 < 0.0f ? 0.0f : y0;
                x1 = x1 > W ? W : x1;
                y1 = y1 > H ? H : y1;
                if (!(x1 > x0 && y1 > y0)) continue;
                const md_gpu_rect_t r = { (uint32_t)x0, (uint32_t)y0, (uint32_t)(x1 - (float)(uint32_t)x0), (uint32_t)(y1 - (float)(uint32_t)y0) };
                if (!r.width || !r.height) continue;
                md_gpu_set_scissor(s, &r);
                a.vertices = va[l] + (md_gpu_addr_t)C->vertex_offset * sizeof(md_post_ui_vertex_t);
                a.tex = C->texture.handle ? C->texture : white;
                a.smp = C->sampler.handle ? C->sampler : p->bilinear;
                ok = md_gpu_draw_indexed(s, p->ui[f][(C->flags & MD_POST_UI_PREMULTIPLIED) ? 1 : 0], ia[l] + (md_gpu_addr_t)C->first_index * isz, d->index_type, C->index_count, 1, &a, sizeof(a));
            }
        }
        ok = md_gpu_render_end(s) && ok;
    }
    if (va) md_free(md_get_heap_allocator(), va, sizeof(md_gpu_addr_t) * 2 * (d->list_count ? d->list_count : 1));
    md_gpu_temp_end(scope);
    if (!ok) MD_LOG_ERROR("md_post_ui: %s", md_gpu_last_error());
    return ok;
}
