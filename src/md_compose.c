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
    md_gpu_addr_t        layers;
    uint32_t             count, stride, width, height;
    uint32_t             _pad0, _pad1;
    md_gpu_storage_tex_t target;
} gpu_compose_args_t;
_Static_assert(sizeof(gpu_compose_args_t) == 192, "Args mirror");

struct md_compose {
    md_gpu_device_t dev;
    md_gpu_kernel_t k_layers;
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
    md_gpu_kernel_desc_t kd = md_shader_md_compose_compose_layers_kernel();
    c->k_layers = md_gpu_kernel_create(device, &kd);
    if (!c->k_layers) {
        MD_LOG_ERROR("md_compose: kernel creation failed: %s", md_gpu_last_error());
        md_compose_destroy(c);
        return NULL;
    }
    return c;
}

void md_compose_destroy(md_compose_ctx_t c) {
    if (!c) return;
    if (c->k_layers) md_gpu_kernel_destroy(c->k_layers);
    md_free(md_get_heap_allocator(), c, sizeof(*c));
}

static md_gpu_float4 f4(float x, float y, float z, float w) {
    md_gpu_float4 v = { x, y, z, w };
    return v;
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
    if (!td || td->type != MD_GPU_TEX_2D || td->format != MD_GPU_FORMAT_RGBA16_FLOAT || !(td->usage & MD_GPU_TEX_STORAGE) ||
        td->width != vd->width || td->height != vd->height) {
        MD_LOG_ERROR("md_compose: the target is to be a 2D STORAGE texture of RGBA16_FLOAT, %u x %u", vd->width, vd->height);
        return false;
    }

    md_float4x4_t P;
    md_memcpy(&P, d->camera.proj, sizeof(P));
    const md_float4x4_t IP = md_float4x4_inverse(P);
    const float* ip = (const float*)&IP;

    const md_compose_shading_t* sh = &d->shading;
    float L[3] = { sh->light_dir[0], sh->light_dir[1], sh->light_dir[2] };
    const float ll = sqrtf(L[0] * L[0] + L[1] * L[1] + L[2] * L[2]);
    if (ll > 0.0f) { L[0] /= ll; L[1] /= ll; L[2] /= ll; } else { L[2] = 1.0f; }
    const bool ortho = d->camera.proj[3] == 0.0f && d->camera.proj[7] == 0.0f && d->camera.proj[11] == 0.0f && d->camera.proj[15] == 1.0f;

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
    a.layers       = t.layers;
    a.count        = t.layers ? t.layer_count : 0;
    a.stride       = t.layer_stride;
    a.width        = vd->width;
    a.height       = vd->height;
    a.target       = md_gpu_texture_storage(d->target, 0);
    md_gpu_grid_t g = md_gpu_grid_for(c->k_layers, vd->width, vd->height, 1);
    return md_gpu_launch(s, c->k_layers, g, &a, sizeof(a));
}
