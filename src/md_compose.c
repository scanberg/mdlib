#include "md_compose.h"

#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_vec_math.h>

#include <math.h>

#include "md_compose_shaders.inl"

/* md_compose.slang's Args */
typedef struct gpu_compose_args_t {
    md_gpu_float4        ip0, ip1, ip2, ip3;
    md_gpu_float4        light, dir_radiance, sky, ground, misc;
    md_gpu_float4        pa, pb, pc;
    md_gpu_addr_t        layers, lin, ao;
    uint32_t             count, stride, width, height;
    uint32_t             level, dir, src, use_ao;
    md_gpu_storage_tex_t target;
    md_gpu_sampled_tex_t depth, color, normal;
    md_gpu_addr_t        vols;
    uint32_t             vol_count, frame;
    md_gpu_sampled_tex_t vf[4];         /* vf0..vf3 */
    md_gpu_sampled_tex_t vt[4];         /* vt0..vt3 */
    md_gpu_sampler_t     trilinear;
    md_gpu_sampled_tex_t vc[4];         /* colour volumes */
    md_gpu_sampled_tex_t vg[4];         /* colour fields */
    md_gpu_sampled_tex_t vm[4];         /* colour maps */
    md_gpu_sampled_tex_t vn[4];         /* normals */
} gpu_compose_args_t;
_Static_assert(sizeof(gpu_compose_args_t) == 496, "Args mirror");

enum { AO_LEVELS = 6, VOL_STRIDE = 36 };   /* md_compose.slang's */

struct md_compose {
    md_gpu_device_t  dev;
    md_gpu_kernel_t  k_linear, k_estimate, k_quality, k_blur, k_compose;
    md_gpu_sampler_t trilinear;
    /* What the compose's volume arguments hold without volumes, never read: a
       null pointer or handle there let lavapipe crash (a load hoisted out of
       a loop that does not run) */
    md_gpu_texture_t none_field, none_transfer, none_color;
};

md_compose_ctx_t md_compose_create(md_gpu_device_t device) {
    if (!device) {
        MD_LOG_ERROR("md_compose: no device");
        return NULL;
    }
    md_compose_ctx_t c = md_alloc(md_get_heap_allocator(), sizeof(*c));
    if (!c) return NULL;
    md_memset(c, 0, sizeof(*c));
    c->dev = device;
    md_gpu_kernel_desc_t kd;
    kd = md_shader_md_compose_ao_linear_kernel();   c->k_linear   = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_compose_ao_estimate_kernel(); c->k_estimate = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_compose_ao_estimate_quality_kernel(); c->k_quality = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_compose_ao_blur_kernel();     c->k_blur     = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_compose_compose_kernel();     c->k_compose  = md_gpu_kernel_create(device, &kd);
    md_gpu_sampler_desc_t sd = {0};
    sd.min_filter = sd.mag_filter = MD_GPU_FILTER_LINEAR;
    c->trilinear = md_gpu_sampler(device, &sd);
    if (!c->k_linear || !c->k_estimate || !c->k_quality || !c->k_blur || !c->k_compose || !c->trilinear.handle) {
        MD_LOG_ERROR("md_compose: kernel creation failed: %s", md_gpu_last_error());
        md_compose_destroy(c);
        return NULL;
    }
    return c;
}

void md_compose_destroy(md_compose_ctx_t c) {
    if (!c) return;
    if (c->k_linear)   md_gpu_kernel_destroy(c->k_linear);
    if (c->k_estimate) md_gpu_kernel_destroy(c->k_estimate);
    if (c->k_quality)  md_gpu_kernel_destroy(c->k_quality);
    if (c->k_blur)     md_gpu_kernel_destroy(c->k_blur);
    if (c->k_compose)  md_gpu_kernel_destroy(c->k_compose);
    if (c->none_field)    md_gpu_texture_destroy(c->none_field);
    if (c->none_transfer) md_gpu_texture_destroy(c->none_transfer);
    if (c->none_color)    md_gpu_texture_destroy(c->none_color);
    md_free(md_get_heap_allocator(), c, sizeof(*c));
}

static md_gpu_float4 f4(float x, float y, float z, float w) {
    md_gpu_float4 v = { x, y, z, w };
    return v;
}

static uint32_t level_dim(uint32_t n, uint32_t l) {
    const uint32_t d = n >> l;
    return d ? d : 1u;
}

/* An optional texture of a volume: 0, or SAMPLED of the type and one of the formats */
static bool vol_tex_ok(md_gpu_texture_t t, md_gpu_tex_type_t type, md_gpu_format_t f0, md_gpu_format_t f1) {
    if (!t) return true;
    const md_gpu_texture_desc_t* d = md_gpu_texture_desc(t);
    return d && d->type == type && (d->usage & MD_GPU_TEX_SAMPLED) && (d->format == f0 || d->format == f1);
}

static bool launch(md_gpu_stream_t s, md_gpu_kernel_t k, uint32_t w, uint32_t h, const gpu_compose_args_t* a) {
    return md_gpu_launch(s, k, md_gpu_grid_for(k, w, h, 1), a, sizeof(*a));
}

bool md_compose(md_gpu_stream_t s, md_compose_ctx_t c, md_render_view_t view, const md_compose_desc_t* d) {
    if (!s || !c || !view || !d) return false;
    const md_render_targets_t t = md_render_view_targets(view);
    const md_gpu_texture_desc_t* vd = t.depth ? md_gpu_texture_desc(t.depth) : NULL;
    const md_gpu_texture_desc_t* td = d->target ? md_gpu_texture_desc(d->target) : NULL;
    if (!vd) {
        MD_LOG_ERROR("md_compose: the view has drawn nothing yet");
        return false;
    }
    const uint32_t W = vd->width, H = vd->height;
    if (!td || td->type != MD_GPU_TEX_2D || td->format != MD_GPU_FORMAT_RGBA16_FLOAT || !(td->usage & MD_GPU_TEX_STORAGE) ||
        td->width != W || td->height != H) {
        MD_LOG_ERROR("md_compose: the target is to be a 2D STORAGE texture of RGBA16_FLOAT, %u x %u", W, H);
        return false;
    }

    const float* P = d->camera.proj;     /* column-major: P[c * 4 + r] */
    md_float4x4_t PM;
    md_memcpy(&PM, P, sizeof(PM));
    const md_float4x4_t IP = md_float4x4_inverse(PM);
    const float* ip = (const float*)&IP;

    const md_compose_shading_t* sh = &d->shading;
    float L[3] = { sh->light_dir[0], sh->light_dir[1], sh->light_dir[2] };
    const float ll = sqrtf(L[0] * L[0] + L[1] * L[1] + L[2] * L[2]);
    if (ll > 0.0f) { L[0] /= ll; L[1] /= ll; L[2] /= ll; } else { L[2] = 1.0f; }
    const bool ortho = P[3] == 0.0f && P[7] == 0.0f && P[11] == 0.0f && P[15] == 1.0f;
    const bool use_ao = d->ao > 0.0f;

    gpu_compose_args_t a;
    md_memset(&a, 0, sizeof(a));
    a.ip0 = f4(ip[0],  ip[1],  ip[2],  ip[3]);
    a.ip1 = f4(ip[4],  ip[5],  ip[6],  ip[7]);
    a.ip2 = f4(ip[8],  ip[9],  ip[10], ip[11]);
    a.ip3 = f4(ip[12], ip[13], ip[14], ip[15]);
    a.light        = f4(L[0], L[1], L[2], sh->light_strength);
    a.dir_radiance = f4(sh->dir_radiance[0], sh->dir_radiance[1], sh->dir_radiance[2], sh->f0);
    a.sky          = f4(sh->sky_radiance[0], sh->sky_radiance[1], sh->sky_radiance[2], sh->roughness);
    a.ground       = f4(sh->ground_radiance[0], sh->ground_radiance[1], sh->ground_radiance[2], sh->ambient);
    a.misc         = f4(sh->specular, sh->mode == MD_COMPOSE_SHADING_BASIC ? 1.0f : 0.0f, sh->linear_colors ? 1.0f : 0.0f, ortho ? 1.0f : 0.0f);
    a.pa           = f4(P[0], P[5], P[8], P[9]);
    a.pb           = f4(P[12], P[13], P[10], P[14]);
    a.pc           = f4(P[11], P[15], 2.0f / (P[5] * (float)H), d->ao);
    a.layers       = t.layers;
    a.count        = t.layers ? t.layer_count : 0;
    a.stride       = t.layer_stride;
    a.width        = W;
    a.height       = H;
    a.use_ao       = use_ao ? 1u : 0u;
    a.target       = md_gpu_texture_storage(d->target, 0);
    a.depth        = md_gpu_texture_sampled(t.depth);
    a.color        = md_gpu_texture_sampled(t.color);
    a.normal       = md_gpu_texture_sampled(t.normal);

    /* The volumes: their records (md_compose.slang's), their textures */
    const uint32_t nv = d->volume_count;
    if (nv > MD_COMPOSE_MAX_VOLUMES || (nv && !d->volumes)) {
        MD_LOG_ERROR("md_compose: at most %d volumes", MD_COMPOSE_MAX_VOLUMES);
        return false;
    }
    md_gpu_float4 rec[MD_COMPOSE_MAX_VOLUMES * VOL_STRIDE];
    md_memset(rec, 0, sizeof(rec));
    for (uint32_t v = 0; v < nv; ++v) {
        const md_compose_volume_t* vol = &d->volumes[v];
        const md_gpu_texture_desc_t* fd = vol->field ? md_gpu_texture_desc(vol->field) : NULL;
        if (!fd || fd->type != MD_GPU_TEX_3D || !(fd->usage & MD_GPU_TEX_SAMPLED) ||
            (fd->format != MD_GPU_FORMAT_R32_FLOAT && fd->format != MD_GPU_FORMAT_R16_FLOAT)) {
            MD_LOG_ERROR("md_compose: volume %u's field is to be a SAMPLED 3D texture of R32_FLOAT or R16_FLOAT", v);
            return false;
        }
        if (!vol_tex_ok(vol->transfer, MD_GPU_TEX_2D, MD_GPU_FORMAT_RGBA8_UNORM, MD_GPU_FORMAT_RGBA8_UNORM) ||
            !vol_tex_ok(vol->color_map, MD_GPU_TEX_2D, MD_GPU_FORMAT_RGBA8_UNORM, MD_GPU_FORMAT_RGBA8_UNORM) ||
            !vol_tex_ok(vol->color_volume, MD_GPU_TEX_3D, MD_GPU_FORMAT_RGBA8_UNORM, MD_GPU_FORMAT_RGBA8_UNORM) ||
            !vol_tex_ok(vol->color_field, MD_GPU_TEX_3D, MD_GPU_FORMAT_R32_FLOAT, MD_GPU_FORMAT_R16_FLOAT) ||
            !vol_tex_ok(vol->normals, MD_GPU_TEX_3D, MD_GPU_FORMAT_RGB10A2_UNORM, MD_GPU_FORMAT_RGBA8_UNORM) ||
            (vol->color_field && !vol->color_map)) {
            MD_LOG_ERROR("md_compose: volume %u: a texture of the wrong type or format, or a colour field without its map (md_compose.h)", v);
            return false;
        }
        if (vol->iso_count > MD_COMPOSE_MAX_ISO || vol->plane_count > MD_COMPOSE_MAX_PLANES || vol->isoline_count > MD_COMPOSE_MAX_ISOLINES) {
            MD_LOG_ERROR("md_compose: volume %u: at most %d isosurfaces, %d planes and %d isolines", v, MD_COMPOSE_MAX_ISO, MD_COMPOSE_MAX_PLANES, MD_COMPOSE_MAX_ISOLINES);
            return false;
        }
        /* view space to its texture space */
        md_float4x4_t V, M;
        md_memcpy(&V, d->camera.view, sizeof(V));
        md_memcpy(&M, vol->model, sizeof(M));
        const md_float4x4_t VT = md_float4x4_inverse(md_float4x4_mul(V, M));
        const float* vt = (const float*)&VT;
        md_gpu_float4* r = rec + v * VOL_STRIDE;
        for (int c4 = 0; c4 < 4; ++c4) r[c4] = f4(vt[4 * c4], vt[4 * c4 + 1], vt[4 * c4 + 2], vt[4 * c4 + 3]);
        const float ext = vol->transfer_range[1] - vol->transfer_range[0];
        r[4] = f4(vol->clip_min[0], vol->clip_min[1], vol->clip_min[2], vol->transfer_range[0]);
        r[5] = f4(vol->clip_max[0], vol->clip_max[1], vol->clip_max[2], ext != 0.0f ? 1.0f / ext : 0.0f);
        r[6] = f4((float)fd->width, (float)fd->height, (float)fd->depth_or_layers, vol->transfer ? 1.0f : 0.0f);
        r[7] = f4((float)vol->iso_count, (float)vol->plane_count, 0.0f, 0.0f);
        float* f = (float*)(r + 8);
        for (uint32_t i = 0; i < vol->iso_count; ++i) {
            f[i]           = vol->iso_values[i];                          /* 8, 9   */
            f[40 + i]      = vol->iso_densities[i];                       /* 18, 19 */
            for (int k = 0; k < 4; ++k) f[8 + 4 * i + k] = vol->iso_colors[i][k];   /* 10..17 */
        }
        for (uint32_t i = 0; i < vol->plane_count; ++i) r[20 + i] = f4(vol->planes[i][0], vol->planes[i][1], vol->planes[i][2], vol->planes[i][3]);
        /* the colour, the isolines, the normals */
        const bool  cvol = vol->color_volume != 0, cfld = !cvol && vol->color_field != 0, nrm = vol->normals != 0;
        const float cext = vol->color_range[1] - vol->color_range[0];
        const uint32_t nl = cfld ? vol->isoline_count : 0u;
        r[28] = f4(vol->color_range[0], cext != 0.0f ? 1.0f / cext : 0.0f, (float)((cvol ? 1u : 0u) | (cfld ? 2u : 0u) | (nrm ? 4u : 0u)),
                   vol->isoline_width > 0.0f ? vol->isoline_width : 1.5f);
        r[29] = f4((float)nl, cfld && vol->isoline_hover ? 1.0f : 0.0f, vol->isoline_hover_value, 0.0f);
        float* il = (float*)(r + 30);
        for (uint32_t i = 0; i < nl; ++i) il[i] = vol->isoline_values[i];
        for (int k = 0; k < 3; ++k) r[32 + k] = f4(d->camera.view[4 * k], d->camera.view[4 * k + 1], d->camera.view[4 * k + 2], 0.0f);
        if (cfld) {
            const md_gpu_texture_desc_t* gd = md_gpu_texture_desc(vol->color_field);
            r[35] = f4((float)gd->width, (float)gd->height, (float)gd->depth_or_layers, 0.0f);
        }
        a.vf[v] = md_gpu_texture_sampled(vol->field);
        if (vol->transfer)     a.vt[v] = md_gpu_texture_sampled(vol->transfer);
        if (cvol)              a.vc[v] = md_gpu_texture_sampled(vol->color_volume);
        if (cfld)              a.vg[v] = md_gpu_texture_sampled(vol->color_field);
        if (cfld)              a.vm[v] = md_gpu_texture_sampled(vol->color_map);
        if (nrm)               a.vn[v] = md_gpu_texture_sampled(vol->normals);
    }
    if (!c->none_field) {
        md_gpu_texture_desc_t nd = {0};
        nd.type = MD_GPU_TEX_3D; nd.format = MD_GPU_FORMAT_R32_FLOAT; nd.usage = MD_GPU_TEX_SAMPLED;
        nd.width = nd.height = nd.depth_or_layers = 1; nd.label = "md_compose.none_field";
        c->none_field = md_gpu_texture_create(s, &nd);
        nd.type = MD_GPU_TEX_3D; nd.format = MD_GPU_FORMAT_RGBA8_UNORM; nd.label = "md_compose.none_color";
        c->none_color = md_gpu_texture_create(s, &nd);
        nd.type = MD_GPU_TEX_2D; nd.format = MD_GPU_FORMAT_RGBA8_UNORM; nd.depth_or_layers = 1; nd.label = "md_compose.none_transfer";
        c->none_transfer = md_gpu_texture_create(s, &nd);
        if (!c->none_field || !c->none_transfer || !c->none_color) {
            MD_LOG_ERROR("md_compose: %s", md_gpu_last_error());
            return false;
        }
    }
    for (uint32_t v = 0; v < 4; ++v) {
        if (!a.vf[v].handle) a.vf[v] = md_gpu_texture_sampled(c->none_field);
        if (!a.vt[v].handle) a.vt[v] = md_gpu_texture_sampled(c->none_transfer);
        if (!a.vc[v].handle) a.vc[v] = md_gpu_texture_sampled(c->none_color);
        if (!a.vg[v].handle) a.vg[v] = md_gpu_texture_sampled(c->none_field);
        if (!a.vm[v].handle) a.vm[v] = md_gpu_texture_sampled(c->none_transfer);
        if (!a.vn[v].handle) a.vn[v] = md_gpu_texture_sampled(c->none_color);
    }
    a.vol_count = nv;
    a.frame     = d->frame;
    a.trilinear = c->trilinear;

    const md_gpu_ordering_t prev = md_gpu_stream_set_ordering(s, MD_GPU_ORDER_IMPLICIT);
    md_gpu_temp_t scope = md_gpu_temp_begin(s);
    bool ok = true;
    a.vols = md_gpu_temp_upload(scope, rec, (size_t)(nv ? nv : 1u) * VOL_STRIDE * sizeof(md_gpu_float4));
    ok = a.vols != 0;
    if (use_ao) {
        /* the linear depth's levels, and two planes at half resolution */
        size_t lin = 0;
        for (uint32_t l = 0; l < AO_LEVELS; ++l) lin += (size_t)level_dim(W, l) * level_dim(H, l);
        const uint32_t hw = level_dim(W, 1), hh = level_dim(H, 1);
        a.lin = md_gpu_temp_alloc(scope, MD_GPU_MEM_DEVICE, lin * sizeof(float)).gpu;
        a.ao  = md_gpu_temp_alloc(scope, MD_GPU_MEM_DEVICE, (size_t)hw * hh * 2u * 8u).gpu;
        ok = a.lin && a.ao;
        for (uint32_t l = 0; ok && l < AO_LEVELS; ++l) {
            a.level = l;
            ok = launch(s, c->k_linear, level_dim(W, l), level_dim(H, l), &a);
        }
        ok = ok && launch(s, d->ao_mode == MD_COMPOSE_AO_PERFORMANCE ? c->k_estimate : c->k_quality, hw, hh, &a);
        a.dir = 0; a.src = 0;
        ok = ok && launch(s, c->k_blur, hw, hh, &a);
        a.dir = 1; a.src = 1;
        ok = ok && launch(s, c->k_blur, hw, hh, &a);   /* back into plane 0 */
    }
    ok = ok && launch(s, c->k_compose, W, H, &a);
    md_gpu_temp_end(scope);
    md_gpu_stream_set_ordering(s, prev);
    if (!ok) MD_LOG_ERROR("md_compose: %s", md_gpu_last_error());
    return ok;
}
