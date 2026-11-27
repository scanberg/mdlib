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
    md_gpu_addr_t        slices;
    uint32_t             slice_count, _sp;
    md_gpu_sampled_tex_t sf[8];         /* the slices' fields */
    md_gpu_sampled_tex_t st[8];         /* and transfer functions */
    md_gpu_float4        lt0, lt1, lt2, lt3;   /* view space to the light's clip space */
    md_gpu_float4        sp;            /* the shadow map's texels across and down, a texel's world size along x, y */
    md_gpu_float4        sq;            /* world per unit of its depth, tan of the light's angular radius, 1: shadows, the contact ray's texels */
    md_gpu_sampled_tex_t shadow;
    float                jit[2];        /* the raster's subpixel offset: the rays through it */
    md_gpu_addr_t        blocks;        /* the volumes' blocks, 1 where one may hold a crossing or emit */
    uint32_t             block_vol, _bp;
    md_gpu_addr_t        tf;            /* the volumes' transfer functions, linear: colour, extinction */
} gpu_compose_args_t;
_Static_assert(sizeof(gpu_compose_args_t) == 784, "Args mirror");

enum { AO_LEVELS = 6, VOL_STRIDE = 37, SLICE_STRIDE = 14, BLOCK = 8 };   /* md_compose.slang's */

struct md_compose {
    md_gpu_device_t  dev;
    md_gpu_kernel_t  k_linear, k_estimate, k_quality, k_blur, k_compose, k_blocks, k_tf;
    md_gpu_sampler_t trilinear;
    /* What the compose's volume arguments hold without volumes, never read: a
       null pointer or handle there let lavapipe crash (a load hoisted out of
       a loop that does not run) */
    md_gpu_texture_t none_field, none_transfer, none_color, none_depth;
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
    kd = md_shader_md_compose_vol_blocks_kernel();  c->k_blocks   = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_compose_vol_transfer_prep_kernel(); c->k_tf = md_gpu_kernel_create(device, &kd);
    md_gpu_sampler_desc_t sd = {0};
    sd.min_filter = sd.mag_filter = MD_GPU_FILTER_LINEAR;
    c->trilinear = md_gpu_sampler(device, &sd);
    if (!c->k_linear || !c->k_estimate || !c->k_quality || !c->k_blur || !c->k_compose || !c->k_blocks || !c->k_tf || !c->trilinear.handle) {
        MD_LOG_ERROR("md_compose: kernel creation failed: %s", md_gpu_last_error());
        md_compose_destroy(c);
        return NULL;
    }
    return c;
}

/* For md_render_bench: the compiled kernels, by index, with what the driver
   reports of each (registers, ...). False past the last. Not in the header. */
bool md_compose_debug_kernel_info(md_compose_ctx_t c, uint32_t i, const char** name, md_gpu_kernel_info_t* info) {
    if (!c) return false;
    const struct { const char* n; md_gpu_kernel_t k; } ks[] = {
        { "compose", c->k_compose }, { "vol_blocks", c->k_blocks }, { "vol_transfer_prep", c->k_tf }, { "ao_linear", c->k_linear }, { "ao_estimate", c->k_estimate },
        { "ao_quality", c->k_quality }, { "ao_blur", c->k_blur },
    };
    if (i >= sizeof(ks) / sizeof(ks[0])) return false;
    if (name) *name = ks[i].n;
    return info ? md_gpu_kernel_info(ks[i].k, info) : true;
}

void md_compose_destroy(md_compose_ctx_t c) {
    if (!c) return;
    if (c->k_linear)   md_gpu_kernel_destroy(c->k_linear);
    if (c->k_estimate) md_gpu_kernel_destroy(c->k_estimate);
    if (c->k_quality)  md_gpu_kernel_destroy(c->k_quality);
    if (c->k_blur)     md_gpu_kernel_destroy(c->k_blur);
    if (c->k_compose)  md_gpu_kernel_destroy(c->k_compose);
    if (c->k_blocks)   md_gpu_kernel_destroy(c->k_blocks);
    if (c->k_tf)       md_gpu_kernel_destroy(c->k_tf);
    if (c->none_field)    md_gpu_texture_destroy(c->none_field);
    if (c->none_transfer) md_gpu_texture_destroy(c->none_transfer);
    if (c->none_color)    md_gpu_texture_destroy(c->none_color);
    if (c->none_depth)    md_gpu_texture_destroy(c->none_depth);
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

static float srgb_decode_c(float x) {
    x = x < 0.0f ? 0.0f : x > 1.0f ? 1.0f : x;
    return x <= 0.04045f ? x / 12.92f : powf((x + 0.055f) / 1.055f, 2.4f);
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
    uint32_t blocks[MD_COMPOSE_MAX_VOLUMES], block_total = 0, texels[MD_COMPOSE_MAX_VOLUMES], texel_total = 0;
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
        const md_gpu_texture_desc_t* td = vol->transfer ? md_gpu_texture_desc(vol->transfer) : NULL;
        texels[v] = td ? td->width : 0;
        r[7] = f4((float)vol->iso_count, (float)vol->plane_count, (float)texels[v], (float)texel_total);
        texel_total += texels[v];
        /* its blocks for the empty space (vol_blocks) */
        const uint32_t gx = (fd->width + BLOCK - 1) / BLOCK, gy = (fd->height + BLOCK - 1) / BLOCK, gz = (fd->depth_or_layers + BLOCK - 1) / BLOCK;
        blocks[v] = gx * gy * gz;
        r[36] = f4((float)gx, (float)gy, (float)gz, (float)block_total);
        block_total += blocks[v];
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
        nd.format = MD_GPU_FORMAT_R32_FLOAT; nd.label = "md_compose.none_depth";
        c->none_depth = md_gpu_texture_create(s, &nd);
        if (!c->none_field || !c->none_transfer || !c->none_color || !c->none_depth) {
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
    /* The slices: their records, their textures */
    const uint32_t ns = d->slice_count;
    if (ns > MD_COMPOSE_MAX_SLICES || (ns && !d->slices)) {
        MD_LOG_ERROR("md_compose: at most %d slices", MD_COMPOSE_MAX_SLICES);
        return false;
    }
    md_gpu_float4 srec[MD_COMPOSE_MAX_SLICES * SLICE_STRIDE];
    md_memset(srec, 0, sizeof(srec));
    for (uint32_t i = 0; i < ns; ++i) {
        const md_compose_slice_t* sl = &d->slices[i];
        /* its edge alone (opacity 0) reads neither */
        const bool shown = sl->opacity > 0.0f;
        if ((shown && (!sl->field || !sl->transfer)) || !vol_tex_ok(sl->field, MD_GPU_TEX_3D, MD_GPU_FORMAT_R32_FLOAT, MD_GPU_FORMAT_R16_FLOAT) ||
            !vol_tex_ok(sl->transfer, MD_GPU_TEX_2D, MD_GPU_FORMAT_RGBA8_UNORM, MD_GPU_FORMAT_RGBA8_UNORM)) {
            MD_LOG_ERROR("md_compose: slice %u: its field is to be a SAMPLED 3D texture of R32_FLOAT or R16_FLOAT, its transfer function a 2D one of RGBA8_UNORM", i);
            return false;
        }
        md_float4x4_t V, R, F;
        md_memcpy(&V, d->camera.view, sizeof(V));
        md_memcpy(&R, sl->rect, sizeof(R));
        md_memcpy(&F, sl->field_from_world, sizeof(F));
        const md_float4x4_t IV = md_float4x4_inverse(V);
        const md_float4x4_t RV = md_float4x4_mul(V, R);                           /* the rectangle's to view */
        const md_float4x4_t VR = md_float4x4_inverse(RV);                         /* view to the rectangle's */
        const md_float4x4_t VF = md_float4x4_mul(F, IV);                          /* view to the field's */
        const float* vr = (const float*)&VR;
        const float* vf = (const float*)&VF;
        const float* rv = (const float*)&RV;
        md_gpu_float4* r = srec + i * SLICE_STRIDE;
        for (int c4 = 0; c4 < 4; ++c4) {
            r[c4]      = f4(vr[4 * c4], vr[4 * c4 + 1], vr[4 * c4 + 2], vr[4 * c4 + 3]);
            r[4 + c4]  = f4(vf[4 * c4], vf[4 * c4 + 1], vf[4 * c4 + 2], vf[4 * c4 + 3]);
            r[10 + c4] = f4(rv[4 * c4], rv[4 * c4 + 1], rv[4 * c4 + 2], rv[4 * c4 + 3]);
        }
        const float ext = sl->transfer_range[1] - sl->transfer_range[0];
        r[8] = f4(sl->transfer_range[0], ext != 0.0f ? 1.0f / ext : 0.0f, sl->opacity, sl->edge_px);
        r[9] = f4(srgb_decode_c(sl->edge[0]), srgb_decode_c(sl->edge[1]), srgb_decode_c(sl->edge[2]), sl->edge[3]);
        a.sf[i] = md_gpu_texture_sampled(sl->field ? sl->field : c->none_field);
        a.st[i] = md_gpu_texture_sampled(sl->transfer ? sl->transfer : c->none_transfer);
    }
    for (uint32_t i = ns; i < MD_COMPOSE_MAX_SLICES; ++i) {
        a.sf[i] = md_gpu_texture_sampled(c->none_field);
        a.st[i] = md_gpu_texture_sampled(c->none_transfer);
    }
    a.slice_count = ns;

    a.vol_count = nv;
    a.frame     = d->frame;
    a.jit[0]    = d->jitter[0];
    a.jit[1]    = d->jitter[1];
    a.trilinear = c->trilinear;

    /* The shadows: the frame's view space to the light's clip space, the map */
    a.shadow = md_gpu_texture_sampled(c->none_depth);
    if (d->shadow) {
        const md_compose_shadow_t* sd = d->shadow;
        const md_gpu_texture_desc_t* md = sd->depth ? md_gpu_texture_desc(sd->depth) : NULL;
        const float* LP = sd->light.proj;
        if (!md || md->type != MD_GPU_TEX_2D || md->format != MD_GPU_FORMAT_D32_FLOAT || !(md->usage & MD_GPU_TEX_SAMPLED) ||
            !(LP[3] == 0.0f && LP[7] == 0.0f && LP[11] == 0.0f && LP[15] == 1.0f) || !(LP[0] > 0.0f && LP[5] > 0.0f && LP[10] > 0.0f)) {
            MD_LOG_ERROR("md_compose: the shadow map is to be a SAMPLED 2D texture of D32_FLOAT, of an orthographic camera (reverse Z)");
            return false;
        }
        md_float4x4_t V, LV, LPM;
        md_memcpy(&V, d->camera.view, sizeof(V));
        md_memcpy(&LV, sd->light.view, sizeof(LV));
        md_memcpy(&LPM, LP, sizeof(LPM));
        const md_float4x4_t VL = md_float4x4_mul(LPM, md_float4x4_mul(LV, md_float4x4_inverse(V)));
        const float* vl = (const float*)&VL;
        a.lt0 = f4(vl[0],  vl[1],  vl[2],  vl[3]);
        a.lt1 = f4(vl[4],  vl[5],  vl[6],  vl[7]);
        a.lt2 = f4(vl[8],  vl[9],  vl[10], vl[11]);
        a.lt3 = f4(vl[12], vl[13], vl[14], vl[15]);
        const float sw = (float)md->width, sh2 = (float)md->height;
        a.sp = f4(sw, sh2, 2.0f / (LP[0] * sw), 2.0f / (LP[5] * sh2));
        a.sq = f4(1.0f / LP[10], tanf(sd->angular_radius > 0.0f ? sd->angular_radius : 0.0f), 1.0f, sd->contact > 0.0f ? sd->contact : 0.0f);
        a.shadow = md_gpu_texture_sampled(sd->depth);
    }

    const md_gpu_ordering_t prev = md_gpu_stream_set_ordering(s, MD_GPU_ORDER_IMPLICIT);
    md_gpu_temp_t scope = md_gpu_temp_begin(s);
    bool ok = true;
    a.vols   = md_gpu_temp_upload(scope, rec, (size_t)(nv ? nv : 1u) * VOL_STRIDE * sizeof(md_gpu_float4));
    a.slices = md_gpu_temp_upload(scope, srec, (size_t)(ns ? ns : 1u) * SLICE_STRIDE * sizeof(md_gpu_float4));
    ok = a.vols != 0 && a.slices != 0;
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
    if (ok && nv) {
        /* the blocks a ray may skip */
        a.blocks = md_gpu_temp_alloc(scope, MD_GPU_MEM_DEVICE, (size_t)block_total * sizeof(uint32_t)).gpu;
        ok = a.blocks != 0;
        for (uint32_t v = 0; ok && v < nv; ++v) {
            a.block_vol = v;
            ok = launch(s, c->k_blocks, blocks[v], 1, &a);
        }
        /* the transfer functions as the march reads them */
        if (ok && texel_total) {
            a.tf = md_gpu_temp_alloc(scope, MD_GPU_MEM_DEVICE, (size_t)texel_total * 16).gpu;
            ok = a.tf != 0;
            for (uint32_t v = 0; ok && v < nv; ++v) {
                if (!texels[v]) continue;
                a.block_vol = v;
                ok = launch(s, c->k_tf, texels[v], 1, &a);
            }
        }
    }
    ok = ok && launch(s, c->k_compose, W, H, &a);
    md_gpu_temp_end(scope);
    md_gpu_stream_set_ordering(s, prev);
    if (!ok) MD_LOG_ERROR("md_compose: %s", md_gpu_last_error());
    return ok;
}
