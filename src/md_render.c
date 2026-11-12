#include "md_render.h"

#include <core/md_allocator.h>
#include <core/md_log.h>

#include <math.h>
#include <stddef.h>
#include <string.h>

#include "md_render_shaders.inl"

/* GPU mirrors of src/shaders/render/md_render.slang. The argument-layout check
   in the shader build holds the Slang side to md_gpu's ABI rule; these are the
   C side of the same rule (md_gpu_float4 at 16-byte offsets, handles at 8). */

enum { M_VIEW, M_PROJ, M_INV_PROJ, M_INV_VIEW, M_VIEW_PREV, M_PROJ_PREV, M_PROJ_CUR, M_COUNT };

typedef struct gpu_view_t {
    md_gpu_float4 m[M_COUNT][4];   /* columns */
    md_gpu_float4 viewport;
    float         dt;
    uint32_t      ortho;
    uint32_t      batch_count;
    uint32_t      total_work;
    float         near_z;
    uint32_t      hiz_valid, hiz_w, hiz_h;
    uint32_t      hiz_levels;
    uint32_t      occlusion, compute_raster;
    float         small_px;
    uint32_t      total_clusters, total_cluster_work, refit_mode, _vpad;
} gpu_view_t;
_Static_assert(sizeof(gpu_view_t) == 528, "ViewData mirror");

enum {
    F_INDEX        = 1u << 0,
    F_VELOCITY     = 1u << 1,
    F_RADIUS       = 1u << 2,
    F_ATTR_INDEX   = 1u << 3,
    F_VISIBLE0     = 1u << 4,
    F_VISIBLE1     = 1u << 5,
    F_STATE0       = 1u << 6,
    F_STATE1       = 1u << 7,
    F_ID           = 1u << 8,
    F_RADIUS_SCALE = 1u << 9,
    F_COLOR        = 1u << 10,
    F_PALETTE_IDX  = 1u << 11,
    F_CLUSTERS     = 1u << 12,
    F_SHARED       = 1u << 13,
};

typedef struct gpu_batch_t {
    uint32_t      kind, count, instance_count, flags;
    uint32_t      stride, vel_stride, attr_count, color_mode;
    md_gpu_addr_t position, velocity, radius, attr_index, color, palette_index, palette;
    md_gpu_addr_t visible0, visible1, state0, state1, index, id, radius_scale;
    md_gpu_addr_t transforms, transforms_prev;
    float         s0, s1;
    uint32_t      uniform_color, palette_count;
    md_gpu_float4 tint;
    float         sharpness, desaturate;
    uint32_t      work_offset, _pad1;
    md_gpu_addr_t cluster_offset, bounds, cstate, _pptr;
    uint32_t      cluster_count, cluster_work_offset, cluster_base, _pad2;
} gpu_batch_t;
_Static_assert(sizeof(gpu_batch_t) == 256, "Batch mirror");

typedef struct gpu_args_t {
    md_gpu_addr_t        view;
    md_gpu_addr_t        batches;
    uint32_t             mode;      /* raster: 0 visibility, 1 mask. cull etc.: phase 1 or 2 */
    uint32_t             list;      /* raster: 0 hw0, 1 hw1, 2 st0, 3 st1 */
    md_gpu_sampled_tex_t id;
    md_gpu_storage_tex_t color, normal, motion, mask;
    uint32_t             width, height, px, py;
    md_gpu_addr_t        pick_out;
    md_gpu_addr_t        lists;
    md_gpu_sampled_tex_t src, hiz;
    md_gpu_storage_tex_t dst;
    uint32_t             src_mip, src_w, src_h, dst_w, dst_h, _pad0;
} gpu_args_t;
_Static_assert(sizeof(gpu_args_t) == 144, "Args mirror");

typedef struct gpu_lists_t {
    md_gpu_addr_t cvis, cocc, cst, hw0, hw1, st0, st1, small, occ, counters, vis64;
} gpu_lists_t;

/* counters[] (md_render.slang C_*) */
enum {
    C_HW1 = 0, C_HW2 = 2, C_ST = 4, C_SMALL1 = 6, C_SMALL2 = 7, C_OCC = 8,
    C_HIDDEN = 9, C_OUTSIDE = 10, C_OCCLUDED = 11,
    C_CVIS1 = 12, C_CVIS2 = 13, C_COCC = 14, C_CST = 15,
    C_CMD = 16, C_STCMD = 32, C_GRID_SMALL = 40, C_GRID_OCC = 48,
    C_GRID_CVIS = 52, C_GRID_COCC = 60, C_GRID_CST = 64, C_CL_VIS = 72, C_CL_OCC = 74,
    C_COUNT = 80,
};
#define IMPLICIT_CLUSTER 64u
#define CTR(v, i) ((v)->counters + (md_gpu_addr_t)(i) * 4u)

#define MAX_BATCHES 4095u
#define PICK_SLOTS  4u
#define STAT_SLOTS  3u
#define CULL_GROUP  64u
#define DEFAULT_SMALL_PIXELS 6.0f

/* =========================================================================
   Context
   ========================================================================= */

struct md_render {
    md_gpu_device_t   dev;
    md_allocator_i*   alloc;
    md_gpu_pipeline_t sphere_vis, sphere_mask, capsule_vis, capsule_mask, merge;
    md_gpu_kernel_t   k_resolve, k_mask_visible, k_pick, k_cull, k_finalize, k_raster_small, k_hiz, k_refit, k_cull_clusters;
    bool              atomic_int64;    /* the compute raster is possible */
};

static md_gpu_pipeline_t make_pipeline(md_gpu_device_t dev, md_gpu_shader_t vs, md_gpu_shader_t fs, bool mask, const char* label) {
    md_gpu_pipeline_desc_t pd = {0};
    pd.vertex       = vs;
    pd.fragment     = fs;
    pd.topology     = MD_GPU_TOPOLOGY_TRIANGLES;
    pd.color_count  = 1;
    pd.depth_format = MD_GPU_FORMAT_D32_FLOAT;
    pd.label        = label;
    if (mask) {
        pd.color[0].format              = MD_GPU_FORMAT_RGBA8_UNORM;
        pd.color[0].blend.enable        = true;
        pd.color[0].blend.src_color     = MD_GPU_BLEND_ONE;
        pd.color[0].blend.dst_color     = MD_GPU_BLEND_ONE;
        pd.color[0].blend.color_op      = MD_GPU_BLEND_OP_MAX;
        pd.color[0].blend.src_alpha     = MD_GPU_BLEND_ONE;
        pd.color[0].blend.dst_alpha     = MD_GPU_BLEND_ONE;
        pd.color[0].blend.alpha_op      = MD_GPU_BLEND_OP_MAX;
        pd.color[0].write_disable       = MD_GPU_COLOR_R | MD_GPU_COLOR_B;   /* the visible part is the kernel's */
    } else {
        pd.color[0].format = MD_GPU_FORMAT_RG32_UINT;
    }
    return md_gpu_pipeline_create(dev, &pd);
}

md_render_t md_render_create(md_gpu_device_t device, const md_render_desc_t* desc) {
    if (!device) {
        MD_LOG_ERROR("md_render: no device");
        return NULL;
    }
    md_gpu_device_info_t info;
    if (!md_gpu_device_info(device, &info) || !info.supports_graphics) {
        MD_LOG_ERROR("md_render: the device has no graphics support");
        return NULL;
    }
    md_allocator_i* alloc = (desc && desc->alloc) ? desc->alloc : md_get_heap_allocator();
    md_render_t r = md_alloc(alloc, sizeof(struct md_render));
    if (!r) return NULL;
    memset(r, 0, sizeof(*r));
    r->dev   = device;
    r->alloc = alloc;

    r->sphere_vis   = make_pipeline(device, md_shader_md_render_sphere_vs_shader(),  md_shader_md_render_sphere_fs_shader(),       false, "md_render.sphere");
    r->sphere_mask  = make_pipeline(device, md_shader_md_render_sphere_vs_shader(),  md_shader_md_render_sphere_mask_fs_shader(),  true,  "md_render.sphere_mask");
    r->capsule_vis  = make_pipeline(device, md_shader_md_render_capsule_vs_shader(), md_shader_md_render_capsule_fs_shader(),      false, "md_render.capsule");
    r->capsule_mask = make_pipeline(device, md_shader_md_render_capsule_vs_shader(), md_shader_md_render_capsule_mask_fs_shader(), true,  "md_render.capsule_mask");
    r->merge        = make_pipeline(device, md_shader_md_render_merge_vs_shader(),   md_shader_md_render_merge_fs_shader(),        false, "md_render.merge");
    r->atomic_int64 = info.supports_atomic_int64;

    md_gpu_kernel_desc_t kd;
    kd = md_shader_md_render_resolve_kernel();      r->k_resolve      = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_mask_visible_kernel(); r->k_mask_visible = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_pick_kernel();         r->k_pick         = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_cull_kernel();         r->k_cull         = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_finalize_kernel();     r->k_finalize     = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_hiz_reduce_kernel();   r->k_hiz          = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_refit_kernel();        r->k_refit        = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_cull_clusters_kernel(); r->k_cull_clusters = md_gpu_kernel_create(device, &kd);
    if (r->atomic_int64) {
        kd = md_shader_md_render_raster_small_kernel(); r->k_raster_small = md_gpu_kernel_create(device, &kd);
        if (!r->k_raster_small) r->atomic_int64 = false;
    }

    if (!r->sphere_vis || !r->sphere_mask || !r->capsule_vis || !r->capsule_mask || !r->merge ||
        !r->k_resolve || !r->k_mask_visible || !r->k_pick || !r->k_cull || !r->k_finalize || !r->k_hiz || !r->k_refit || !r->k_cull_clusters) {
        MD_LOG_ERROR("md_render: pipeline or kernel creation failed: %s", md_gpu_last_error());
        md_render_destroy(r);
        return NULL;
    }
    return r;
}

void md_render_destroy(md_render_t r) {
    if (!r) return;
    if (r->sphere_vis)     md_gpu_pipeline_destroy(r->sphere_vis);
    if (r->sphere_mask)    md_gpu_pipeline_destroy(r->sphere_mask);
    if (r->capsule_vis)    md_gpu_pipeline_destroy(r->capsule_vis);
    if (r->capsule_mask)   md_gpu_pipeline_destroy(r->capsule_mask);
    if (r->k_resolve)      md_gpu_kernel_destroy(r->k_resolve);
    if (r->k_mask_visible) md_gpu_kernel_destroy(r->k_mask_visible);
    if (r->k_pick)         md_gpu_kernel_destroy(r->k_pick);
    if (r->merge)          md_gpu_pipeline_destroy(r->merge);
    if (r->k_cull)         md_gpu_kernel_destroy(r->k_cull);
    if (r->k_finalize)     md_gpu_kernel_destroy(r->k_finalize);
    if (r->k_raster_small) md_gpu_kernel_destroy(r->k_raster_small);
    if (r->k_hiz)          md_gpu_kernel_destroy(r->k_hiz);
    if (r->k_refit)        md_gpu_kernel_destroy(r->k_refit);
    if (r->k_cull_clusters) md_gpu_kernel_destroy(r->k_cull_clusters);
    md_free(r->alloc, r, sizeof(struct md_render));
}

/* =========================================================================
   Small math (double where it is cheap and matters: inverses)
   ========================================================================= */

static void mat_mul(float out[16], const float a[16], const float b[16]) {
    float t[16];
    for (int c = 0; c < 4; ++c)
        for (int r = 0; r < 4; ++r) {
            double s = 0;
            for (int k = 0; k < 4; ++k) s += (double)a[k * 4 + r] * (double)b[c * 4 + k];
            t[c * 4 + r] = (float)s;
        }
    memcpy(out, t, sizeof(t));
}

static bool mat_inverse(float out[16], const float m[16]) {
    double a[16], inv[16];
    for (int i = 0; i < 16; ++i) a[i] = m[i];
    inv[0]  =  a[5]*a[10]*a[15] - a[5]*a[11]*a[14] - a[9]*a[6]*a[15] + a[9]*a[7]*a[14] + a[13]*a[6]*a[11] - a[13]*a[7]*a[10];
    inv[4]  = -a[4]*a[10]*a[15] + a[4]*a[11]*a[14] + a[8]*a[6]*a[15] - a[8]*a[7]*a[14] - a[12]*a[6]*a[11] + a[12]*a[7]*a[10];
    inv[8]  =  a[4]*a[9]*a[15]  - a[4]*a[11]*a[13] - a[8]*a[5]*a[15] + a[8]*a[7]*a[13] + a[12]*a[5]*a[11] - a[12]*a[7]*a[9];
    inv[12] = -a[4]*a[9]*a[14]  + a[4]*a[10]*a[13] + a[8]*a[5]*a[14] - a[8]*a[6]*a[13] - a[12]*a[5]*a[10] + a[12]*a[6]*a[9];
    inv[1]  = -a[1]*a[10]*a[15] + a[1]*a[11]*a[14] + a[9]*a[2]*a[15] - a[9]*a[3]*a[14] - a[13]*a[2]*a[11] + a[13]*a[3]*a[10];
    inv[5]  =  a[0]*a[10]*a[15] - a[0]*a[11]*a[14] - a[8]*a[2]*a[15] + a[8]*a[3]*a[14] + a[12]*a[2]*a[11] - a[12]*a[3]*a[10];
    inv[9]  = -a[0]*a[9]*a[15]  + a[0]*a[11]*a[13] + a[8]*a[1]*a[15] - a[8]*a[3]*a[13] - a[12]*a[1]*a[11] + a[12]*a[3]*a[9];
    inv[13] =  a[0]*a[9]*a[14]  - a[0]*a[10]*a[13] - a[8]*a[1]*a[14] + a[8]*a[2]*a[13] + a[12]*a[1]*a[10] - a[12]*a[2]*a[9];
    inv[2]  =  a[1]*a[6]*a[15]  - a[1]*a[7]*a[14]  - a[5]*a[2]*a[15] + a[5]*a[3]*a[14] + a[13]*a[2]*a[7]  - a[13]*a[3]*a[6];
    inv[6]  = -a[0]*a[6]*a[15]  + a[0]*a[7]*a[14]  + a[4]*a[2]*a[15] - a[4]*a[3]*a[14] - a[12]*a[2]*a[7]  + a[12]*a[3]*a[6];
    inv[10] =  a[0]*a[5]*a[15]  - a[0]*a[7]*a[13]  - a[4]*a[1]*a[15] + a[4]*a[3]*a[13] + a[12]*a[1]*a[7]  - a[12]*a[3]*a[5];
    inv[14] = -a[0]*a[5]*a[14]  + a[0]*a[6]*a[13]  + a[4]*a[1]*a[14] - a[4]*a[2]*a[13] - a[12]*a[1]*a[6]  + a[12]*a[2]*a[5];
    inv[3]  = -a[1]*a[6]*a[11]  + a[1]*a[7]*a[10]  + a[5]*a[2]*a[11] - a[5]*a[3]*a[10] - a[9]*a[2]*a[7]   + a[9]*a[3]*a[6];
    inv[7]  =  a[0]*a[6]*a[11]  - a[0]*a[7]*a[10]  - a[4]*a[2]*a[11] + a[4]*a[3]*a[10] + a[8]*a[2]*a[7]   - a[8]*a[3]*a[6];
    inv[11] = -a[0]*a[5]*a[11]  + a[0]*a[7]*a[9]   + a[4]*a[1]*a[11] - a[4]*a[3]*a[9]  - a[8]*a[1]*a[7]   + a[8]*a[3]*a[5];
    inv[15] =  a[0]*a[5]*a[10]  - a[0]*a[6]*a[9]   - a[4]*a[1]*a[10] + a[4]*a[2]*a[9]  + a[8]*a[1]*a[6]   - a[8]*a[2]*a[5];
    double det = a[0] * inv[0] + a[1] * inv[4] + a[2] * inv[8] + a[3] * inv[12];
    if (det == 0.0 || !isfinite(det)) return false;
    det = 1.0 / det;
    for (int i = 0; i < 16; ++i) out[i] = (float)(inv[i] * det);
    return true;
}

static void put_mat(md_gpu_float4 dst[4], const float m[16]) {
    for (int c = 0; c < 4; ++c) {
        dst[c].x = m[c * 4 + 0];
        dst[c].y = m[c * 4 + 1];
        dst[c].z = m[c * 4 + 2];
        dst[c].w = m[c * 4 + 3];
    }
}

static const float IDENTITY[16] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };

static float halton(uint32_t i, uint32_t base) {
    float f = 1.0f, r = 0.0f;
    while (i > 0) {
        f /= (float)base;
        r += f * (float)(i % base);
        i /= base;
    }
    return r;
}

/* =========================================================================
   Views
   ========================================================================= */

typedef struct pick_meta_t {
    uint64_t key;
    uint32_t kind;
    uint32_t count;
} pick_meta_t;

typedef struct pick_slot_t {
    bool          pending;
    uint64_t      seq;
    uint32_t      x, y;
    md_gpu_addr_t dev;             /* device scratch the kernel writes */
    md_gpu_mem_t  host;            /* HOST_READ copy                    */
    md_gpu_sync_t sync;
    pick_meta_t*  meta;            /* the batches of the raster it read */
    uint32_t      meta_count;
} pick_slot_t;

struct md_render_view {
    md_render_t         r;
    char                label[64];

    uint32_t            width, height;
    md_render_targets_t t;

    /* The frame of the last draw, and of the last raster */
    bool                has_raster;
    md_render_frame_t   last;              /* scalars; batches below           */
    md_render_batch_t*  batches;           /* [batch_cap]                      */
    float*              xforms;            /* transforms of `batches`, packed  */
    uint32_t*           xform_offset;      /* [batch_cap + 1] in floats        */
    uint32_t            batch_cap;
    uint32_t            xform_cap;         /* floats                           */
    pick_meta_t*        raster_meta;       /* batches as of the last raster    */
    uint32_t            raster_meta_count;

    double              raster_time;
    float               raster_view[16];
    float               raster_proj[16];   /* unjittered                       */
    uint32_t            sample;            /* of the targets                   */
    float               jitter[2];

    /* GPU */
    md_gpu_addr_t       view_buf;
    md_gpu_addr_t       batch_buf;
    uint32_t            batch_buf_cap;
    md_gpu_addr_t       xform_buf;         /* current then previous per batch  */
    size_t              xform_buf_cap;     /* bytes                            */
    gpu_batch_t*        table;             /* host copy of the batch table     */

    /* Cull lists, counters, compute raster buffer, depth pyramid */
    md_gpu_addr_t       lists_buf;         /* gpu_lists_t                      */
    md_gpu_addr_t       list_mem;          /* the lists themselves             */
    size_t              list_mem_cap;
    md_gpu_addr_t       counters;          /* uint32 [C_COUNT]                 */
    md_gpu_addr_t       vis64;
    size_t              vis64_cap;
    md_gpu_addr_t       bounds_mem;        /* per cluster: 2 float4 bounds, 1 uint state */
    size_t              bounds_cap;
    bool                bounds_valid;      /* of the last raster's geometry   */
    bool                state_lists_ok;    /* hold the current state bits     */
    md_gpu_texture_t    hiz;
    uint32_t            hiz_w, hiz_h, hiz_levels;
    bool                hiz_ok;            /* holds the last raster            */
    uint32_t            raster_count;

    struct {
        md_gpu_mem_t    host;
        md_gpu_sync_t   sync;
        bool            pending;
        uint32_t        raster;
        uint64_t        elements;
        uint32_t        clusters;
        bool            refitted;
    } stats[STAT_SLOTS];
    uint32_t            stats_returned;

    pick_slot_t         picks[PICK_SLOTS];
    uint64_t            pick_seq;
    uint64_t            pick_returned;
};

md_render_view_t md_render_view_create(md_render_t r, const char* label) {
    if (!r) return NULL;
    md_render_view_t v = md_alloc(r->alloc, sizeof(struct md_render_view));
    if (!v) return NULL;
    memset(v, 0, sizeof(*v));
    v->r = r;
    if (label) {
        strncpy(v->label, label, sizeof(v->label) - 1);
    }
    return v;
}

static void destroy_targets(md_render_view_t v) {
    md_gpu_texture_t* ts[] = { &v->t.depth, &v->t.id, &v->t.color, &v->t.normal, &v->t.motion, &v->t.mask };
    for (size_t i = 0; i < sizeof(ts) / sizeof(ts[0]); ++i) {
        if (*ts[i]) md_gpu_texture_destroy(*ts[i]);
        *ts[i] = NULL;
    }
}

static void free_pick(md_render_view_t v, md_gpu_stream_t s, pick_slot_t* p) {
    if (p->dev)      md_gpu_free(s, p->dev);
    if (p->host.gpu) md_gpu_free(s, p->host.gpu);
    if (p->meta)     md_free(v->r->alloc, p->meta, p->meta_count * sizeof(pick_meta_t));
    memset(p, 0, sizeof(*p));
}

void md_render_view_destroy(md_gpu_stream_t s, md_render_view_t v) {
    if (!v) return;
    md_allocator_i* a = v->r->alloc;
    destroy_targets(v);
    if (v->view_buf)  md_gpu_free(s, v->view_buf);
    if (v->batch_buf) md_gpu_free(s, v->batch_buf);
    if (v->xform_buf) md_gpu_free(s, v->xform_buf);
    if (v->lists_buf) md_gpu_free(s, v->lists_buf);
    if (v->bounds_mem) md_gpu_free(s, v->bounds_mem);
    if (v->list_mem)  md_gpu_free(s, v->list_mem);
    if (v->counters)  md_gpu_free(s, v->counters);
    if (v->vis64)     md_gpu_free(s, v->vis64);
    if (v->hiz)       md_gpu_texture_destroy(v->hiz);
    for (uint32_t i = 0; i < STAT_SLOTS; ++i) if (v->stats[i].host.gpu) md_gpu_free(s, v->stats[i].host.gpu);
    for (uint32_t i = 0; i < PICK_SLOTS; ++i) free_pick(v, s, &v->picks[i]);
    if (v->batches)      md_free(a, v->batches, v->batch_cap * sizeof(md_render_batch_t));
    if (v->table)        md_free(a, v->table, v->batch_cap * sizeof(gpu_batch_t));
    if (v->xform_offset) md_free(a, v->xform_offset, (v->batch_cap + 1) * sizeof(uint32_t));
    if (v->xforms)       md_free(a, v->xforms, v->xform_cap * sizeof(float));
    if (v->raster_meta)  md_free(a, v->raster_meta, v->raster_meta_count * sizeof(pick_meta_t));
    md_free(a, v, sizeof(struct md_render_view));
}

static md_gpu_texture_t make_target(md_gpu_stream_t s, md_gpu_format_t fmt, md_gpu_tex_usage_t usage, uint32_t w, uint32_t h, const char* label) {
    md_gpu_texture_desc_t td = {0};
    td.type   = MD_GPU_TEX_2D;
    td.format = fmt;
    td.usage  = usage;
    td.width  = w;
    td.height = h;
    td.label  = label;
    return md_gpu_texture_create(s, &td);
}

static bool ensure_targets(md_gpu_stream_t s, md_render_view_t v, uint32_t w, uint32_t h) {
    if (v->t.depth && v->width == w && v->height == h) return true;
    destroy_targets(v);
    const md_gpu_tex_usage_t RT = MD_GPU_TEX_RENDER_TARGET, SA = MD_GPU_TEX_SAMPLED, ST = MD_GPU_TEX_STORAGE;
    v->t.depth  = make_target(s, MD_GPU_FORMAT_D32_FLOAT,   RT | SA,      w, h, "md_render.depth");
    v->t.id     = make_target(s, MD_GPU_FORMAT_RG32_UINT,   RT | SA | ST, w, h, "md_render.id");
    v->t.color  = make_target(s, MD_GPU_FORMAT_RGBA8_UNORM, SA | ST,      w, h, "md_render.color");
    v->t.normal = make_target(s, MD_GPU_FORMAT_RG16_FLOAT,  SA | ST,      w, h, "md_render.normal");
    v->t.motion = make_target(s, MD_GPU_FORMAT_RG16_FLOAT,  SA | ST,      w, h, "md_render.motion");
    v->t.mask   = make_target(s, MD_GPU_FORMAT_RGBA8_UNORM, RT | SA | ST, w, h, "md_render.mask");
    if (!v->t.depth || !v->t.id || !v->t.color || !v->t.normal || !v->t.motion || !v->t.mask) {
        MD_LOG_ERROR("md_render: target creation failed: %s", md_gpu_last_error());
        destroy_targets(v);
        return false;
    }
    v->width  = w;
    v->height = h;
    return true;
}

/* ---- Frame comparison ------------------------------------------------------- */

/* The part of a batch that rasterisation depends on, with padding and the
   transform pointer zeroed so that memcmp compares meaning. */
static void canonical_geometry(md_render_batch_t* out, const md_render_batch_t* b) {
    memset(out, 0, sizeof(*out));
    out->key    = b->key;
    out->kind   = b->kind;
    out->count  = b->count;
    out->points.count           = b->points.count;
    out->points.position        = b->points.position;
    out->points.stride          = b->points.stride;
    out->points.velocity_stride = b->points.velocity_stride;
    out->points.velocity        = b->points.velocity;
    out->points.radius          = b->points.radius;
    out->points.attr_count      = b->points.attr_count;
    out->points.attr_index      = b->points.attr_index;
    out->points.visible[0]      = b->points.visible[0];
    out->points.visible[1]      = b->points.visible[1];
    out->points.state[0]        = b->points.state[0];   /* addresses: the mask pass draws from them */
    out->points.state[1]        = b->points.state[1];
    out->index           = b->index;
    out->id              = b->id;
    out->cluster_offset  = b->cluster_offset;
    out->cluster_count   = b->cluster_offset ? b->cluster_count : 0;
    if (b->kind == MD_RENDER_SPHERES) {
        out->style.spheres.radius_scale = b->style.spheres.radius_scale;
        out->style.spheres.radius_add   = b->style.spheres.radius_add;
    } else if (b->kind == MD_RENDER_CAPSULES) {
        const md_render_capsule_style_t* c = &b->style.capsules;
        out->style.capsules.radius       = c->radius;
        out->style.capsules.cap[0]       = c->cap[0];
        out->style.capsules.cap[1]       = c->cap[1];
        out->style.capsules.radius_scale = c->radius_scale;
        out->style.capsules.dash_length  = c->dash_length;
        out->style.capsules.gap_length   = c->gap_length;
        out->style.capsules.offset       = c->offset;
        out->style.capsules.offset_mode  = c->offset_mode;
        memcpy(out->style.capsules.cell, c->cell, sizeof(c->cell));
    }
    out->transform_count = b->transform ? b->transform_count : 0;
}

static void canonical_color(md_render_color_t* out, const md_render_batch_t* b) {
    memset(out, 0, sizeof(*out));
    out->mode          = b->color.mode;
    out->uniform       = b->color.uniform;
    out->palette       = b->color.palette;
    out->palette_count = b->color.palette_count;
    out->sharpness     = b->color.sharpness;
    memcpy(out->tint, b->color.tint, sizeof(out->tint));
    out->desaturate    = b->color.desaturate;
}

static bool validate(const md_render_frame_t* f) {
    if (f->width == 0 || f->height == 0) {
        MD_LOG_ERROR("md_render: frame of zero size");
        return false;
    }
    if (f->batch_count > MAX_BATCHES) {
        MD_LOG_ERROR("md_render: %u batches; at most %u per frame", f->batch_count, MAX_BATCHES);
        return false;
    }
    if (f->batch_count && !f->batches) {
        MD_LOG_ERROR("md_render: batch_count without batches");
        return false;
    }
    uint64_t total = 0;
    for (uint32_t i = 0; i < f->batch_count; ++i) total += (f->batches[i].kind == MD_RENDER_SPHERES || f->batches[i].kind == MD_RENDER_CAPSULES)
        ? (uint64_t)f->batches[i].count * ((f->batches[i].transform && f->batches[i].transform_count) ? f->batches[i].transform_count : 1u) : 0u;
    if (total >= (1ull << 31)) {
        MD_LOG_ERROR("md_render: %llu elements times instances in a frame; at most 2^31", (unsigned long long)total);
        return false;
    }
    for (uint32_t i = 0; i < f->batch_count; ++i) {
        const md_render_batch_t* b = &f->batches[i];
        if (b->key == 0) {
            MD_LOG_ERROR("md_render: batch %u has key 0", i);
            return false;
        }
        for (uint32_t j = 0; j < i; ++j) {
            if (f->batches[j].key == b->key) {
                MD_LOG_ERROR("md_render: batches %u and %u share the key %llu", j, i, (unsigned long long)b->key);
                return false;
            }
        }
        switch (b->kind) {
        case MD_RENDER_SPHERES:
        case MD_RENDER_CAPSULES:
            break;
        case MD_RENDER_EXTERNAL:
            continue;
        case MD_RENDER_POLYGONS:
        case MD_RENDER_SPLINES:
        case MD_RENDER_TRIANGLES:
            MD_LOG_ERROR("md_render: batch %u: kind %d is not implemented yet", i, (int)b->kind);
            return false;
        default:
            MD_LOG_ERROR("md_render: batch %u: invalid kind %d", i, (int)b->kind);
            return false;
        }
        if (b->count == 0) continue;
        if (!b->points.position) {
            MD_LOG_ERROR("md_render: batch %u has no positions", i);
            return false;
        }
        if (b->cluster_offset && b->cluster_count == 0) {
            MD_LOG_ERROR("md_render: batch %u has cluster offsets but no cluster count", i);
            return false;
        }
        if ((b->points.stride & 3u) || (b->points.velocity_stride & 3u)) {
            MD_LOG_ERROR("md_render: batch %u: strides must be multiples of 4", i);
            return false;
        }
        if (b->kind == MD_RENDER_CAPSULES) {
            const md_render_capsule_style_t* c = &b->style.capsules;
            if (c->cap[0] != MD_RENDER_CAP_ROUND || c->cap[1] != MD_RENDER_CAP_ROUND || c->dash_length != 0.0f || c->offset) {
                MD_LOG_ERROR("md_render: batch %u: flat caps, dashes and offsets are not implemented yet", i);
                return false;
            }
            for (int k = 0; k < 9; ++k) {
                if (c->cell[k] != 0.0f) {
                    MD_LOG_ERROR("md_render: batch %u: periodic capsules are not implemented yet", i);
                    return false;
                }
            }
        }
        uint64_t inst = b->transform ? b->transform_count : 1;
        if (inst == 0) inst = 1;
        if (inst * (uint64_t)b->count >= (1ull << 50)) {
            MD_LOG_ERROR("md_render: batch %u: more than 2^50 elements times instances", i);
            return false;
        }
    }
    return true;
}

static uint32_t instance_count(const md_render_batch_t* b) {
    return (b->transform && b->transform_count) ? b->transform_count : 1u;
}

/* Work items of the cull: elements times instances, for the kinds that are culled */
static uint64_t batch_work(const md_render_batch_t* b) {
    if (b->kind != MD_RENDER_SPHERES && b->kind != MD_RENDER_CAPSULES) return 0;
    return (uint64_t)b->count * instance_count(b);
}

/* Store the frame as the one last drawn. */
/* Clusters of a batch: given, or runs of 64 consecutive elements */
static uint32_t batch_clusters(const md_render_batch_t* b) {
    if (b->kind != MD_RENDER_SPHERES && b->kind != MD_RENDER_CAPSULES) return 0;
    if (b->cluster_offset && b->cluster_count) return b->cluster_count;
    return (b->count + IMPLICIT_CLUSTER - 1) / IMPLICIT_CLUSTER;
}

/* Two batches whose cluster bounds are the same refit: the same points,
   index, clusters, visibility, state and style; not the same key, ids,
   velocities, colours or instances */
static bool bounds_shared(const md_render_batch_t* a, const md_render_batch_t* b) {
    md_render_batch_t x, y;
    canonical_geometry(&x, a);
    canonical_geometry(&y, b);
    x.key = y.key = 0;
    x.id = y.id = 0;
    x.transform_count = y.transform_count = 0;
    x.points.velocity = y.points.velocity = 0;
    x.points.velocity_stride = y.points.velocity_stride = 0;
    return memcmp(&x, &y, sizeof(x)) == 0;
}

static bool ensure_batch_cap(md_render_view_t v, uint32_t n) {
    md_allocator_i* a = v->r->alloc;
    if (n <= v->batch_cap && v->batches) return true;
    uint32_t cap = n < 16 ? 16 : n;
    v->batches      = md_realloc(a, v->batches,      v->batch_cap * sizeof(md_render_batch_t), cap * sizeof(md_render_batch_t));
    v->table        = md_realloc(a, v->table,        v->batch_cap * sizeof(gpu_batch_t),       cap * sizeof(gpu_batch_t));
    v->xform_offset = md_realloc(a, v->xform_offset, (v->batch_cap ? v->batch_cap + 1 : 0) * sizeof(uint32_t), (cap + 1) * sizeof(uint32_t));
    if (!v->batches || !v->table || !v->xform_offset) return false;
    if (v->batch_cap == 0) v->xform_offset[0] = 0;
    v->batch_cap = cap;
    return true;
}

static bool remember(md_render_view_t v, const md_render_frame_t* f) {
    md_allocator_i* a = v->r->alloc;
    if (!ensure_batch_cap(v, f->batch_count)) return false;
    uint32_t floats = 0;
    for (uint32_t i = 0; i < f->batch_count; ++i) floats += 16u * instance_count(&f->batches[i]);
    if (floats > v->xform_cap) {
        v->xforms = md_realloc(a, v->xforms, v->xform_cap * sizeof(float), floats * sizeof(float));
        if (!v->xforms) return false;
        v->xform_cap = floats;
    }
    uint32_t o = 0;
    for (uint32_t i = 0; i < f->batch_count; ++i) {
        const md_render_batch_t* b = &f->batches[i];
        v->batches[i] = *b;
        v->batches[i].transform = NULL;
        v->xform_offset[i] = o;
        uint32_t n = instance_count(b);
        if (b->transform && b->transform_count) memcpy(v->xforms + o, b->transform, n * 16 * sizeof(float));
        else memcpy(v->xforms + o, IDENTITY, sizeof(IDENTITY));
        o += 16u * n;
    }
    v->xform_offset[f->batch_count] = o;
    v->last = *f;
    v->last.batches = NULL;
    return true;
}

typedef enum change_t { CHANGE_NONE = 0, CHANGE_STATE, CHANGE_APPEARANCE, CHANGE_RASTER } change_t;

static change_t compare(const md_render_view_t v, const md_render_frame_t* f) {
    if (!v->has_raster) return CHANGE_RASTER;
    const md_render_frame_t* l = &v->last;
    if (f->width != l->width || f->height != l->height) return CHANGE_RASTER;
    if (memcmp(&f->camera, &l->camera, sizeof(f->camera)) != 0) return CHANGE_RASTER;
    if (f->batch_count != l->batch_count) return CHANGE_RASTER;
    if (f->content[MD_RENDER_CONTENT_GEOMETRY] != l->content[MD_RENDER_CONTENT_GEOMETRY]) return CHANGE_RASTER;
    if (f->options != l->options || f->small_pixels != l->small_pixels) return CHANGE_RASTER;
    change_t c = CHANGE_NONE;
    for (uint32_t i = 0; i < f->batch_count; ++i) {
        md_render_batch_t x, y;
        canonical_geometry(&x, &f->batches[i]);
        canonical_geometry(&y, &v->batches[i]);
        if (memcmp(&x, &y, sizeof(x)) != 0) return CHANGE_RASTER;
        uint32_t n = instance_count(&f->batches[i]);
        const float* t = (f->batches[i].transform && f->batches[i].transform_count) ? f->batches[i].transform : IDENTITY;
        if (memcmp(t, v->xforms + v->xform_offset[i], n * 16 * sizeof(float)) != 0) return CHANGE_RASTER;
        md_render_color_t cx, cy;
        canonical_color(&cx, &f->batches[i]);
        canonical_color(&cy, &v->batches[i]);
        if (memcmp(&cx, &cy, sizeof(cx)) != 0) c = CHANGE_APPEARANCE;
    }
    if (f->content[MD_RENDER_CONTENT_APPEARANCE] != l->content[MD_RENDER_CONTENT_APPEARANCE]) c = CHANGE_APPEARANCE;
    if (c == CHANGE_NONE && f->content[MD_RENDER_CONTENT_STATE] != l->content[MD_RENDER_CONTENT_STATE]) c = CHANGE_STATE;
    return c;
}

/* The cluster bounds of the last raster hold for this frame: same batches
   (but for colours and instance transforms) and no geometry announced */
static bool geometry_unchanged(const md_render_view_t v, const md_render_frame_t* f) {
    if (!v->has_raster || !v->bounds_valid) return false;
    if (f->batch_count != v->last.batch_count) return false;
    if (f->content[MD_RENDER_CONTENT_GEOMETRY] != v->last.content[MD_RENDER_CONTENT_GEOMETRY]) return false;
    for (uint32_t i = 0; i < f->batch_count; ++i) {
        md_render_batch_t x, y;
        canonical_geometry(&x, &f->batches[i]);
        canonical_geometry(&y, &v->batches[i]);
        x.transform_count = y.transform_count = 0;
        if (memcmp(&x, &y, sizeof(x)) != 0) return false;
    }
    return true;
}

/* ---- GPU tables --------------------------------------------------------------- */

static md_gpu_addr_t align_addr(md_gpu_addr_t x) { return (x + 255u) & ~(md_gpu_addr_t)255u; }

static bool gpu_realloc(md_gpu_stream_t s, md_gpu_addr_t* addr, size_t* cap, size_t need) {
    if (*addr && *cap >= need) return true;
    if (*addr) md_gpu_free(s, *addr);
    size_t c = need < 4096 ? 4096 : need + need / 2;
    *addr = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, c).gpu;
    *cap  = *addr ? c : 0;
    return *addr != 0;
}

static void fill_batch(gpu_batch_t* g, const md_render_batch_t* b, uint32_t n_inst, md_gpu_addr_t xf, md_gpu_addr_t xf_prev) {
    memset(g, 0, sizeof(*g));
    const md_render_points_t* p = &b->points;
    g->kind           = (uint32_t)b->kind;
    g->count          = b->count;
    g->instance_count = n_inst;
    g->stride         = p->stride ? p->stride : 12u;
    g->vel_stride     = p->velocity_stride ? p->velocity_stride : 8u;
    g->attr_count     = p->attr_count ? p->attr_count : p->count;
    g->color_mode     = (uint32_t)b->color.mode;
    g->position       = p->position;
    g->velocity       = p->velocity;
    g->radius         = p->radius;
    g->attr_index     = p->attr_index;
    g->color          = p->color;
    g->palette_index  = p->palette_index;
    g->palette        = b->color.palette;
    g->visible0       = p->visible[0];
    g->visible1       = p->visible[1];
    g->state0         = p->state[0];
    g->state1         = p->state[1];
    g->index          = b->index;
    g->id             = b->id;
    g->transforms      = xf;
    g->transforms_prev = xf_prev;
    uint32_t f = 0;
    if (b->index)          f |= F_INDEX;
    if (p->velocity)       f |= F_VELOCITY;
    if (p->radius)         f |= F_RADIUS;
    if (p->attr_index)     f |= F_ATTR_INDEX;
    if (p->visible[0])     f |= F_VISIBLE0;
    if (p->visible[1])     f |= F_VISIBLE1;
    if (p->state[0])       f |= F_STATE0;
    if (p->state[1])       f |= F_STATE1;
    if (b->id)             f |= F_ID;
    if (p->color)          f |= F_COLOR;
    if (p->palette_index)  f |= F_PALETTE_IDX;
    if (b->kind == MD_RENDER_SPHERES) {
        g->s0 = b->style.spheres.radius_scale;
        g->s1 = b->style.spheres.radius_add;
    } else if (b->kind == MD_RENDER_CAPSULES) {
        g->s0 = b->style.capsules.radius;
        g->radius_scale = b->style.capsules.radius_scale;
        if (b->style.capsules.radius_scale) f |= F_RADIUS_SCALE;
    }
    g->flags         = f;
    g->uniform_color = b->color.uniform;
    g->palette_count = b->color.palette_count;
    g->tint.x = b->color.tint[0];
    g->tint.y = b->color.tint[1];
    g->tint.z = b->color.tint[2];
    g->tint.w = b->color.tint[3];
    g->sharpness  = b->color.sharpness;
    g->desaturate = b->color.desaturate;
}

/* The previous transforms of batch i: those of the batch with its key at the
   last raster, when it had as many instances; otherwise its current ones. */
static const float* prev_transforms(const md_render_view_t v, const md_render_frame_t* prev_frame_scalars,
                                    const md_render_batch_t* prev_batches, const float* prev_xforms,
                                    const uint32_t* prev_offset, const md_render_batch_t* b, const float* cur) {
    if (!prev_frame_scalars) return cur;
    (void)v;
    for (uint32_t j = 0; j < prev_frame_scalars->batch_count; ++j) {
        if (prev_batches[j].key == b->key) {
            uint32_t n = (prev_offset[j + 1] - prev_offset[j]) / 16u;
            return n == instance_count(b) ? prev_xforms + prev_offset[j] : cur;
        }
    }
    return cur;
}

/* Builds and uploads the batch table and transforms for `f`. At a raster the
   previous transforms come from the frame last rasterised; otherwise the
   transform buffer is left as it is (geometry is unchanged). */
static bool upload_tables(md_gpu_stream_t s, md_render_view_t v, const md_render_frame_t* f, bool raster) {
    if (raster) {
        uint32_t floats = 0;
        for (uint32_t i = 0; i < f->batch_count; ++i) floats += 16u * instance_count(&f->batches[i]);
        size_t bytes = (size_t)floats * 2u * sizeof(float);
        if (bytes == 0) bytes = 64;
        float* host = md_alloc(v->r->alloc, bytes);
        if (!host) return false;
        /* Previous transforms are read from what remember() kept of the last raster */
        const md_render_frame_t* pf = v->has_raster ? &v->last : NULL;
        uint32_t o = 0;
        for (uint32_t i = 0; i < f->batch_count; ++i) {
            const md_render_batch_t* b = &f->batches[i];
            uint32_t n = instance_count(b);
            const float* cur = (b->transform && b->transform_count) ? b->transform : IDENTITY;
            const float* prv = prev_transforms(v, pf, v->batches, v->xforms, v->xform_offset, b, cur);
            memcpy(host + o,               cur, n * 16 * sizeof(float));
            memcpy(host + floats + o,      prv, n * 16 * sizeof(float));
            o += 16u * n;
        }
        if (!gpu_realloc(s, &v->xform_buf, &v->xform_buf_cap, bytes)) {
            md_free(v->r->alloc, host, bytes);
            return false;
        }
        md_gpu_upload(s, v->xform_buf, host, bytes);
        md_free(v->r->alloc, host, bytes);
        /* table addresses */
        o = 0;
        for (uint32_t i = 0; i < f->batch_count; ++i) {
            uint32_t n = instance_count(&f->batches[i]);
            md_gpu_addr_t xf = v->xform_buf + (md_gpu_addr_t)o * sizeof(float);
            md_gpu_addr_t xp = v->xform_buf + (md_gpu_addr_t)(floats + o) * sizeof(float);
            v->table[i].transforms      = xf;   /* filled again below, kept for resolves */
            v->table[i].transforms_prev = xp;
            o += 16u * n;
        }
    }
    /* Cluster bounds: 2 float4 per cluster, then a uint of state bits each */
    uint32_t clusters = 0;
    for (uint32_t i = 0; i < f->batch_count; ++i) clusters += batch_clusters(&f->batches[i]);
    {
        size_t need = (size_t)(clusters ? clusters : 1) * 40u + 256u;
        md_gpu_addr_t old = v->bounds_mem;
        if (!gpu_realloc(s, &v->bounds_mem, &v->bounds_cap, need)) return false;
        if (v->bounds_mem != old) v->bounds_valid = false;
    }
    const md_gpu_addr_t cstate_base = align_addr(v->bounds_mem + (md_gpu_addr_t)clusters * 32u);   /* 2 uint per cluster */
    uint32_t work = 0, cbase = 0, cwork = 0;
    for (uint32_t i = 0; i < f->batch_count; ++i) {
        const md_render_batch_t* b = &f->batches[i];
        md_gpu_addr_t xf = v->table[i].transforms, xp = v->table[i].transforms_prev;
        fill_batch(&v->table[i], b, instance_count(b), xf, xp);
        gpu_batch_t* g = &v->table[i];
        g->work_offset         = work;
        g->cluster_count       = batch_clusters(b);
        g->cluster_offset      = (b->cluster_offset && b->cluster_count) ? b->cluster_offset : 0;
        if (g->cluster_offset) g->flags |= F_CLUSTERS;
        g->cluster_base        = cbase;
        g->cluster_work_offset = cwork;
        g->bounds              = v->bounds_mem + (md_gpu_addr_t)cbase * 32u;
        g->cstate              = cstate_base + (md_gpu_addr_t)cbase * 8u;
        for (uint32_t j = 0; j < i; ++j) {
            if (g->cluster_count && bounds_shared(b, &f->batches[j])) {
                g->bounds = v->table[j].bounds;
                g->cstate = v->table[j].cstate;
                g->flags |= F_SHARED;
                break;
            }
        }
        work  += (uint32_t)batch_work(b);
        cbase += g->cluster_count;
        cwork += g->cluster_count * instance_count(b);
    }
    size_t cap = (size_t)v->batch_buf_cap * sizeof(gpu_batch_t);
    size_t need = (size_t)(f->batch_count ? f->batch_count : 1) * sizeof(gpu_batch_t);
    if (!gpu_realloc(s, &v->batch_buf, &cap, need)) return false;
    v->batch_buf_cap = (uint32_t)(cap / sizeof(gpu_batch_t));
    if (f->batch_count) md_gpu_upload(s, v->batch_buf, v->table, f->batch_count * sizeof(gpu_batch_t));
    return true;
}

static gpu_args_t base_args(const md_render_view_t v) {
    gpu_args_t a;
    memset(&a, 0, sizeof(a));
    a.view    = v->view_buf;
    a.batches = v->batch_buf;
    a.id      = md_gpu_texture_sampled(v->t.id);
    a.color   = md_gpu_texture_storage(v->t.color, 0);
    a.normal  = md_gpu_texture_storage(v->t.normal, 0);
    a.motion  = md_gpu_texture_storage(v->t.motion, 0);
    a.mask    = md_gpu_texture_storage(v->t.mask, 0);
    a.width   = v->width;
    a.height  = v->height;
    a.lists   = v->lists_buf;
    if (v->hiz) a.hiz = md_gpu_texture_sampled(v->hiz);
    return a;
}

/* ---- Cull lists, pyramid -------------------------------------------------------- */


/* Lists sized for the worst case of this frame: every element in each list of
   its kind. 32 bytes per element in all; clusters would cut that down. */
static bool ensure_lists(md_gpu_stream_t s, md_render_view_t v, uint64_t work0, uint64_t work1, uint64_t cwork, bool compute) {
    const uint64_t tot = work0 + work1;
    const uint64_t sizes[9] = { cwork, cwork, cwork, work0, work1, work0, work1, tot, tot };
    size_t need = 0;
    for (int i = 0; i < 9; ++i) need += (size_t)align_addr((sizes[i] + 1) * 8u);
    if (!gpu_realloc(s, &v->list_mem, &v->list_mem_cap, need)) return false;
    gpu_lists_t l;
    md_gpu_addr_t o = v->list_mem;
    md_gpu_addr_t* dst[9] = { &l.cvis, &l.cocc, &l.cst, &l.hw0, &l.hw1, &l.st0, &l.st1, &l.small, &l.occ };
    for (int i = 0; i < 9; ++i) { *dst[i] = o; o += align_addr((sizes[i] + 1) * 8u); }
    if (!v->counters) v->counters = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, C_COUNT * 4u).gpu;
    if (!v->lists_buf) v->lists_buf = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, sizeof(gpu_lists_t)).gpu;
    if (!v->counters || !v->lists_buf) return false;
    if (compute && !gpu_realloc(s, &v->vis64, &v->vis64_cap, (size_t)v->width * v->height * 8u)) return false;
    l.counters = v->counters;
    l.vis64    = v->vis64;
    return md_gpu_upload(s, v->lists_buf, &l, sizeof(l));
}

/* Level 0 is half the target, rounded up, down to 1 x 1 */
static bool ensure_hiz(md_gpu_stream_t s, md_render_view_t v) {
    uint32_t w = (v->width + 1) / 2, h = (v->height + 1) / 2;
    if (v->hiz && v->hiz_w == w && v->hiz_h == h) return true;
    if (v->hiz) md_gpu_texture_destroy(v->hiz);
    uint32_t levels = 1, m = w > h ? w : h;
    while (m > 1) { m >>= 1; ++levels; }
    md_gpu_texture_desc_t td = {0};
    td.type       = MD_GPU_TEX_2D;
    td.format     = MD_GPU_FORMAT_R32_FLOAT;
    td.usage      = MD_GPU_TEX_SAMPLED | MD_GPU_TEX_STORAGE;
    td.width      = w;
    td.height     = h;
    td.mip_levels = levels;
    td.label      = "md_render.hiz";
    v->hiz        = md_gpu_texture_create(s, &td);
    v->hiz_w      = w;
    v->hiz_h      = h;
    v->hiz_levels = levels;
    v->hiz_ok     = false;
    return v->hiz != NULL;
}

static bool build_hiz(md_gpu_stream_t s, md_render_view_t v) {
    gpu_args_t a = base_args(v);
    uint32_t sw = v->width, sh = v->height, dw = v->hiz_w, dh = v->hiz_h;
    for (uint32_t L = 0; L < v->hiz_levels; ++L) {
        a.src     = L == 0 ? md_gpu_texture_sampled(v->t.depth) : md_gpu_texture_sampled(v->hiz);
        a.src_mip = L == 0 ? 0 : L - 1;
        a.src_w   = sw;
        a.src_h   = sh;
        a.dst     = md_gpu_texture_storage(v->hiz, L);
        a.dst_w   = dw;
        a.dst_h   = dh;
        if (!md_gpu_launch(s, v->r->k_hiz, md_gpu_grid_for(v->r->k_hiz, dw, dh, 1), &a, sizeof(a))) return false;
        sw = dw; sh = dh;
        dw = dw > 1 ? (dw + 1) / 2 : 1;
        dh = dh > 1 ? (dh + 1) / 2 : 1;
    }
    return true;
}

/* One pass of a cull kernel: over n items directly, or over a grid a
   finalize wrote (n == 0, grid at counter word `grid`) */
static bool pass(md_gpu_stream_t s, md_render_view_t v, md_gpu_kernel_t k, uint32_t mode, uint32_t n, uint32_t grid) {
    gpu_args_t a = base_args(v);
    a.mode = mode;
    if (grid) return md_gpu_launch_indirect(s, k, CTR(v, grid), &a, sizeof(a));
    uint32_t groups = (n + CULL_GROUP - 1) / CULL_GROUP;
    if (groups > 65535) groups = 65535;
    if (groups == 0) groups = 1;
    return md_gpu_launch(s, k, md_gpu_grid(groups, 1, 1), &a, sizeof(a));
}

static bool fin(md_gpu_stream_t s, md_render_view_t v, uint32_t mode) {
    gpu_args_t a = base_args(v);
    a.mode = mode;
    return md_gpu_launch(s, v->r->k_finalize, md_gpu_grid(1, 1, 1), &a, sizeof(a));
}

/* The state lists, for the masks: clusters with a state bit in the frustum,
   then their elements with one, occluded or not */
static bool run_state_lists(md_gpu_stream_t s, md_render_view_t v, uint32_t cluster_work) {
    md_render_t r = v->r;
    return pass(s, v, r->k_cull_clusters, 3, cluster_work, 0) && fin(s, v, 31) &&
           pass(s, v, r->k_cull, 3, 0, C_GRID_CST) && fin(s, v, 32);
}

static bool run_raster_small(md_gpu_stream_t s, md_render_view_t v, uint32_t phase) {
    gpu_args_t a = base_args(v);
    a.mode = phase;
    return md_gpu_launch_indirect(s, v->r->k_raster_small, CTR(v, C_GRID_SMALL + (phase - 1) * 4), &a, sizeof(a));
}

/* Phase 1 clears; phase 2 adds. The compute raster's result goes in first, as
   depth and ids, then each kind in one indirect draw whose instance count the
   cull wrote. The CPU's cost is independent of the number of batches. */
static bool run_raster(md_gpu_stream_t s, md_render_view_t v, uint32_t phase, bool compute, const bool has_kind[2]) {
    md_render_t r = v->r;
    md_gpu_render_desc_t rd = {0};
    rd.color_count       = 1;
    rd.color[0].texture  = v->t.id;
    rd.color[0].load     = phase == 1 ? MD_GPU_LOAD_CLEAR : MD_GPU_LOAD;
    rd.color[0].store    = MD_GPU_STORE;
    rd.depth.texture     = v->t.depth;
    rd.depth.load        = phase == 1 ? MD_GPU_LOAD_CLEAR : MD_GPU_LOAD;
    rd.depth.store       = MD_GPU_STORE;
    rd.depth.clear_depth = 0.0f;
    rd.label             = phase == 1 ? "md_render.raster1" : "md_render.raster2";
    if (!md_gpu_render_begin(s, &rd)) return false;
    gpu_args_t a = base_args(v);
    md_gpu_draw_state_t ds = {0};
    ds.depth_compare = MD_GPU_COMPARE_GREATER_EQUAL;
    ds.depth_write   = true;
    if (compute) {
        md_gpu_set_draw_state(s, &ds);
        md_gpu_draw(s, r->merge, 3, 1, &a, sizeof(a));
    }
    for (uint32_t k = 0; k < 2; ++k) {
        if (!has_kind[k]) continue;
        ds.cull = k == 0 ? MD_GPU_CULL_NONE : MD_GPU_CULL_BACK;
        md_gpu_set_draw_state(s, &ds);
        a.mode = 0;
        a.list = k;
        md_gpu_draw_indirect(s, k == 0 ? r->sphere_vis : r->capsule_vis, CTR(v, C_CMD + ((phase - 1) * 2 + k) * 4), 1, &a, sizeof(a));
    }
    return md_gpu_render_end(s);
}

static bool run_resolve(md_gpu_stream_t s, md_render_view_t v) {
    gpu_args_t a = base_args(v);
    md_gpu_grid_t g = md_gpu_grid_for(v->r->k_resolve, v->width, v->height, 1);
    return md_gpu_launch(s, v->r->k_resolve, g, &a, sizeof(a));
}

/* The visible part from the ids; the hidden part by rasterising the state
   lists of the cull without a depth test, so it costs as much as the
   selection and not the system. */
static bool run_masks(md_gpu_stream_t s, md_render_view_t v, const md_render_frame_t* f) {
    gpu_args_t a = base_args(v);
    md_gpu_grid_t g = md_gpu_grid_for(v->r->k_mask_visible, v->width, v->height, 1);
    if (!md_gpu_launch(s, v->r->k_mask_visible, g, &a, sizeof(a))) return false;

    bool has[2] = { false, false };
    for (uint32_t i = 0; i < f->batch_count; ++i) {
        const md_render_batch_t* b = &f->batches[i];
        if (b->count && (b->points.state[0] || b->points.state[1])) {
            if (b->kind == MD_RENDER_SPHERES)  has[0] = true;
            if (b->kind == MD_RENDER_CAPSULES) has[1] = true;
        }
    }
    if (!has[0] && !has[1]) return true;
    md_gpu_render_desc_t rd = {0};
    rd.color_count      = 1;
    rd.color[0].texture = v->t.mask;
    rd.color[0].load    = MD_GPU_LOAD;
    rd.color[0].store   = MD_GPU_STORE;
    rd.depth.texture    = v->t.depth;     /* not tested: the fragment stage writes depth, which Metal wants an attachment for */
    rd.depth.load       = MD_GPU_LOAD;
    rd.depth.store      = MD_GPU_STORE;
    rd.label            = "md_render.mask";
    if (!md_gpu_render_begin(s, &rd)) return false;
    md_gpu_draw_state_t ds = {0};
    ds.depth_compare = MD_GPU_COMPARE_ALWAYS;
    for (uint32_t k = 0; k < 2; ++k) {
        if (!has[k]) continue;
        ds.cull = k == 0 ? MD_GPU_CULL_NONE : MD_GPU_CULL_BACK;
        md_gpu_set_draw_state(s, &ds);
        a.mode = 1;
        a.list = 2 + k;
        md_gpu_draw_indirect(s, k == 0 ? v->r->sphere_mask : v->r->capsule_mask, CTR(v, C_STCMD + k * 4), 1, &a, sizeof(a));
    }
    return md_gpu_render_end(s);
}

static void record_stats(md_gpu_stream_t s, md_render_view_t v, uint64_t elements, uint32_t clusters, bool refitted) {
    /* the slot not pending, else the oldest */
    uint32_t best = 0;
    for (uint32_t i = 0; i < STAT_SLOTS; ++i) {
        if (!v->stats[i].pending) { best = i; break; }
        if (v->stats[i].raster < v->stats[best].raster) best = i;
    }
    if (!v->stats[best].host.gpu) v->stats[best].host = md_gpu_malloc(s, MD_GPU_MEM_HOST_READ, C_COUNT * 4u);
    if (!v->stats[best].host.cpu) return;
    md_gpu_copy(s, v->stats[best].host.gpu, v->counters, C_COUNT * 4u);
    v->stats[best].sync     = md_gpu_stream_record(s);
    v->stats[best].pending  = true;
    v->stats[best].raster   = v->raster_count;
    v->stats[best].elements = elements;
    v->stats[best].clusters = clusters;
    v->stats[best].refitted = refitted;
}

bool md_render_view_stats(md_render_view_t v, md_render_stats_t* out) {
    if (!v || !out) return false;
    int best = -1;
    for (uint32_t i = 0; i < STAT_SLOTS; ++i) {
        if (!v->stats[i].pending || !md_gpu_sync_is_complete(v->stats[i].sync)) continue;
        if (v->stats[i].raster <= v->stats_returned) { v->stats[i].pending = false; continue; }
        if (best < 0 || v->stats[i].raster > v->stats[best].raster) best = (int)i;
    }
    if (best < 0) return false;
    const uint32_t* c = (const uint32_t*)v->stats[best].host.cpu;
    memset(out, 0, sizeof(*out));
    out->elements    = v->stats[best].elements;
    out->hidden      = c[C_HIDDEN];
    out->outside     = c[C_OUTSIDE];
    out->occluded    = c[C_OCCLUDED];
    out->hardware[0] = c[C_HW1] + c[C_HW1 + 1];
    out->hardware[1] = c[C_HW2] + c[C_HW2 + 1];
    out->compute[0]  = c[C_SMALL1];
    out->compute[1]  = c[C_SMALL2];
    out->occlusion_retested = c[C_OCC];
    out->clusters             = v->stats[best].clusters;
    out->clusters_expanded[0] = c[C_CVIS1];
    out->clusters_expanded[1] = c[C_CVIS2];
    out->clusters_occluded    = c[C_CL_OCC];
    out->refitted             = v->stats[best].refitted;
    out->raster      = v->stats[best].raster;
    v->stats_returned = out->raster;
    for (uint32_t i = 0; i < STAT_SLOTS; ++i) if (v->stats[i].raster <= v->stats_returned) v->stats[i].pending = false;
    return true;
}

static void store_raster_meta(md_render_view_t v, const md_render_frame_t* f) {
    md_allocator_i* a = v->r->alloc;
    if (v->raster_meta) md_free(a, v->raster_meta, v->raster_meta_count * sizeof(pick_meta_t));
    v->raster_meta = NULL;
    v->raster_meta_count = f->batch_count;
    if (!f->batch_count) return;
    v->raster_meta = md_alloc(a, f->batch_count * sizeof(pick_meta_t));
    if (!v->raster_meta) { v->raster_meta_count = 0; return; }
    for (uint32_t i = 0; i < f->batch_count; ++i) {
        v->raster_meta[i].key   = f->batches[i].key;
        v->raster_meta[i].kind  = (uint32_t)f->batches[i].kind;
        v->raster_meta[i].count = f->batches[i].count;
    }
}

bool md_render_view_draw(md_gpu_stream_t s, md_render_view_t v, const md_render_frame_t* f, md_render_status_t* out) {
    if (out) memset(out, 0, sizeof(*out));
    if (!s || !v || !f || !validate(f)) return false;

    if (!ensure_batch_cap(v, f->batch_count)) return false;
    bool resized = !(v->t.depth && v->width == f->width && v->height == f->height);
    if (!ensure_targets(s, v, f->width, f->height)) return false;

    change_t c = resized ? CHANGE_RASTER : compare(v, f);
    const bool geometry_same = geometry_unchanged(v, f);
    const bool state_changed = v->has_raster && f->content[MD_RENDER_CONTENT_STATE] != v->last.content[MD_RENDER_CONTENT_STATE];
    uint32_t samples = f->samples ? f->samples : 1u;

    bool raster = false, resolve = false, mask = false;
    uint32_t sample = v->sample;
    switch (c) {
    case CHANGE_RASTER:
        raster = resolve = mask = true;
        sample = 0;
        break;
    case CHANGE_APPEARANCE:
        resolve = mask = true;
        sample = 0;
        break;
    case CHANGE_STATE:
        mask = true;
        break;
    case CHANGE_NONE:
        if (f->jitter == MD_RENDER_JITTER_ALWAYS ||
            (f->jitter == MD_RENDER_JITTER_STATIC && v->sample + 1 < samples)) {
            raster = resolve = mask = true;
            sample = v->sample + 1;
        }
        break;
    }

    const bool occlusion = !(f->options & MD_RENDER_NO_OCCLUSION_CULLING);
    const bool compute   = !(f->options & MD_RENDER_NO_COMPUTE_RASTER) && v->r->atomic_int64;
    bool hiz_valid = false;
    if (resized) v->hiz_ok = false;

    if (raster) {
        /* Jitter, in pixels. Sample 0 is pixel centres, so a view drawn once
           (or a test) sees exactly what its camera says. */
        float jx = 0.0f, jy = 0.0f;
        if (f->jitter != MD_RENDER_JITTER_NONE && sample > 0) {
            jx = halton(sample, 2) - 0.5f;
            jy = halton(sample, 3) - 0.5f;
        }
        float proj[16];
        memcpy(proj, f->camera.proj, sizeof(proj));
        const float ndx = 2.0f * jx / (float)f->width, ndy = -2.0f * jy / (float)f->height;
        for (int col = 0; col < 4; ++col) {
            proj[col * 4 + 0] += ndx * proj[col * 4 + 3];
            proj[col * 4 + 1] += ndy * proj[col * 4 + 3];
        }
        gpu_view_t gv;
        memset(&gv, 0, sizeof(gv));
        float inv_proj[16], inv_view[16];
        if (!mat_inverse(inv_proj, proj) || !mat_inverse(inv_view, f->camera.view)) {
            MD_LOG_ERROR("md_render: singular camera matrix");
            return false;
        }
        put_mat(gv.m[M_VIEW],      f->camera.view);
        put_mat(gv.m[M_PROJ],      proj);
        put_mat(gv.m[M_INV_PROJ],  inv_proj);
        put_mat(gv.m[M_INV_VIEW],  inv_view);
        put_mat(gv.m[M_VIEW_PREV], v->has_raster ? v->raster_view : f->camera.view);
        put_mat(gv.m[M_PROJ_PREV], v->has_raster ? v->raster_proj : f->camera.proj);
        put_mat(gv.m[M_PROJ_CUR],  f->camera.proj);
        gv.viewport.x = (float)f->width;
        gv.viewport.y = (float)f->height;
        gv.viewport.z = 1.0f / (float)f->width;
        gv.viewport.w = 1.0f / (float)f->height;
        gv.dt          = v->has_raster ? (float)(f->time - v->raster_time) : 0.0f;
        gv.ortho       = (f->camera.proj[3] == 0.0f && f->camera.proj[7] == 0.0f && f->camera.proj[11] == 0.0f && f->camera.proj[15] == 1.0f) ? 1u : 0u;
        gv.batch_count = f->batch_count;
        {
            uint64_t tw = 0;
            for (uint32_t i = 0; i < f->batch_count; ++i) tw += batch_work(&f->batches[i]);
            gv.total_work = (uint32_t)tw;
        }
        /* The near plane: where the projection's depth is 1 */
        {
            const float* P = f->camera.proj;
            float den = P[10] - P[11];
            gv.near_z = den != 0.0f ? (P[15] - P[14]) / den : 0.0f;
        }
        if (occlusion && !ensure_hiz(s, v)) return false;
        hiz_valid          = occlusion && v->hiz_ok;
        gv.hiz_valid       = hiz_valid ? 1u : 0u;
        gv.hiz_w           = v->hiz_w;
        gv.hiz_h           = v->hiz_h;
        gv.hiz_levels      = v->hiz_levels;
        gv.occlusion       = occlusion ? 1u : 0u;
        gv.compute_raster  = compute ? 1u : 0u;
        gv.small_px        = f->small_pixels > 0.0f ? f->small_pixels : DEFAULT_SMALL_PIXELS;
        if (!v->view_buf) v->view_buf = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, sizeof(gpu_view_t)).gpu;
        if (!v->view_buf) return false;
        md_gpu_upload(s, v->view_buf, &gv, sizeof(gv));

        if (!upload_tables(s, v, f, true)) return false;
        const bool refit = !geometry_same;
        gv.total_clusters     = 0;
        gv.total_cluster_work = 0;
        for (uint32_t i = 0; i < f->batch_count; ++i) {
            gv.total_clusters     += batch_clusters(&f->batches[i]);
            gv.total_cluster_work += batch_clusters(&f->batches[i]) * instance_count(&f->batches[i]);
        }
        gv.refit_mode = refit ? 0u : 1u;
        md_gpu_upload(s, v->view_buf, &gv, sizeof(gv));

        uint64_t work[2] = { 0, 0 };
        bool has_kind[2] = { false, false };
        for (uint32_t i = 0; i < f->batch_count; ++i) {
            const md_render_batch_t* b = &f->batches[i];
            uint32_t k = b->kind == MD_RENDER_CAPSULES ? 1u : 0u;
            work[k] += batch_work(b);
            if (batch_work(b)) has_kind[k] = true;
        }
        const uint32_t total = (uint32_t)(work[0] + work[1]);
        const md_gpu_addr_t list_mem_before = v->list_mem;
        if (!ensure_lists(s, v, work[0], work[1], gv.total_cluster_work, compute)) return false;
        /* State lists: rebuilt only when the state bits or the geometry
           changed (they do not depend on the camera); their counts and
           commands survive the counters' reset otherwise */
        const bool state_lists = refit || state_changed || !v->state_lists_ok || v->list_mem != list_mem_before;
        if (state_lists) {
            md_gpu_memset(s, v->counters, 0, C_COUNT * 4u);
        } else {
            static const uint32_t keep[][2] = { { C_ST, C_ST + 2 }, { C_CST, C_CST + 1 }, { C_STCMD, C_STCMD + 8 }, { C_GRID_CST, C_GRID_CST + 3 } };
            uint32_t w = 0;
            for (size_t i = 0; i < sizeof(keep) / sizeof(keep[0]); ++i) {
                if (keep[i][0] > w) md_gpu_memset(s, CTR(v, w), 0, (keep[i][0] - w) * 4u);
                w = keep[i][1];
            }
            if (w < C_COUNT) md_gpu_memset(s, CTR(v, w), 0, (C_COUNT - w) * 4u);
        }
        if (compute) md_gpu_memset(s, v->vis64, 0, (size_t)v->width * v->height * 8u);

        /* Cluster bounds: refitted when the geometry changed, else only their
           state bits when those did. Not when only the camera moved. */
        md_render_t r = v->r;
        if ((refit || state_changed) && !pass(s, v, r->k_refit, 0, gv.total_clusters, 0)) return false;
        v->bounds_valid = true;

        /* Phase 1: clusters, then the elements of the visible ones, against
           the previous raster's pyramid (when there is one) */
        if (!pass(s, v, r->k_cull_clusters, 1, gv.total_cluster_work, 0) || !fin(s, v, 11)) return false;
        if (!pass(s, v, r->k_cull, 1, 0, C_GRID_CVIS) || !fin(s, v, 12)) return false;
        if (compute && !run_raster_small(s, v, 1)) return false;
        if (!run_raster(s, v, 1, compute, has_kind)) return false;

        if (occlusion) {
            if (!build_hiz(s, v)) return false;
            if (hiz_valid) {
                /* Phase 2: the clusters and the elements phase 1 found
                   occluded, against this raster */
                if (!pass(s, v, r->k_cull_clusters, 2, 0, C_GRID_COCC) || !fin(s, v, 21)) return false;
                if (!pass(s, v, r->k_cull, 2, 0, C_GRID_OCC)) return false;
                if (!pass(s, v, r->k_cull, 4, 0, C_GRID_CVIS + 4) || !fin(s, v, 22)) return false;
                if (compute && !run_raster_small(s, v, 2)) return false;
                if (!run_raster(s, v, 2, compute, has_kind)) return false;
                if (!build_hiz(s, v)) return false;     /* for the next raster */
            }
            v->hiz_ok = true;
        } else {
            v->hiz_ok = false;
        }
        if (state_lists) {
            if (!run_state_lists(s, v, gv.total_cluster_work)) return false;
            v->state_lists_ok = true;
        }
        v->raster_count++;
        record_stats(s, v, total, gv.total_cluster_work, refit);

        v->raster_time = f->time;
        memcpy(v->raster_view, f->camera.view, sizeof(v->raster_view));
        memcpy(v->raster_proj, f->camera.proj, sizeof(v->raster_proj));
        v->jitter[0] = jx;
        v->jitter[1] = jy;
        store_raster_meta(v, f);
    } else if (resolve) {
        if (!upload_tables(s, v, f, false)) return false;
    }
    if (!raster && state_changed) {
        /* Only the state bits changed: refit the clusters' state bits, then
           the state lists from them. Cost: the clusters, and the elements of
           the clusters with a state bit. */
        uint32_t clusters = 0, cwork = 0;
        for (uint32_t i = 0; i < f->batch_count; ++i) {
            clusters += batch_clusters(&f->batches[i]);
            cwork    += batch_clusters(&f->batches[i]) * instance_count(&f->batches[i]);
        }
        const uint32_t one = 1;
        md_gpu_upload(s, v->view_buf + offsetof(gpu_view_t, refit_mode), &one, 4);
        if (!pass(s, v, v->r->k_refit, 0, clusters, 0)) return false;
        md_gpu_memset(s, CTR(v, C_ST), 0, 8u);
        md_gpu_memset(s, CTR(v, C_CST), 0, 4u);
        if (!run_state_lists(s, v, cwork)) return false;
    }

    if (resolve && !run_resolve(s, v)) return false;
    if (mask && !run_masks(s, v, f)) return false;

    if (!remember(v, f)) return false;
    v->has_raster = true;
    v->sample     = sample;

    if (out) {
        out->work = (raster ? MD_RENDER_WORK_RASTER : 0u) | (resolve ? MD_RENDER_WORK_RESOLVE : 0u) | (mask ? MD_RENDER_WORK_MASK : 0u);
        out->sample    = sample;
        out->jitter[0] = v->jitter[0];
        out->jitter[1] = v->jitter[1];
    }
    return true;
}

md_render_targets_t md_render_view_targets(md_render_view_t v) {
    md_render_targets_t t = {0};
    if (v) t = v->t;
    return t;
}

void md_render_encode_id(uint32_t batch_index, uint32_t part, uint64_t linear, uint32_t out_id[2]) {
    out_id[0] = (uint32_t)(linear & 0xFFFFFFFFu);
    out_id[1] = ((batch_index + 1u) << 20) | ((part & 3u) << 18) | (uint32_t)((linear >> 32) & 0x3FFFFu);
}

/* ---- Picking -------------------------------------------------------------------- */

bool md_render_view_pick(md_gpu_stream_t s, md_render_view_t v, uint32_t x, uint32_t y) {
    if (!s || !v || !v->has_raster) return false;
    /* A free slot, else the oldest */
    pick_slot_t* p = NULL;
    for (uint32_t i = 0; i < PICK_SLOTS; ++i) {
        if (!v->picks[i].pending) { p = &v->picks[i]; break; }
    }
    if (!p) {
        p = &v->picks[0];
        for (uint32_t i = 1; i < PICK_SLOTS; ++i) if (v->picks[i].seq < p->seq) p = &v->picks[i];
    }
    if (!p->dev)      p->dev  = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, 64).gpu;
    if (!p->host.gpu) p->host = md_gpu_malloc(s, MD_GPU_MEM_HOST_READ, 64);
    if (!p->dev || !p->host.cpu) return false;

    md_allocator_i* a = v->r->alloc;
    if (p->meta) md_free(a, p->meta, p->meta_count * sizeof(pick_meta_t));
    p->meta = NULL;
    p->meta_count = v->raster_meta_count;
    if (p->meta_count) {
        p->meta = md_alloc(a, p->meta_count * sizeof(pick_meta_t));
        if (!p->meta) { p->meta_count = 0; return false; }
        memcpy(p->meta, v->raster_meta, p->meta_count * sizeof(pick_meta_t));
    }

    gpu_args_t args = base_args(v);
    args.px = x;
    args.py = y;
    args.pick_out = p->dev;
    md_gpu_grid_t g = { 1, 1, 1 };
    if (!md_gpu_launch(s, v->r->k_pick, g, &args, sizeof(args))) return false;
    if (!md_gpu_copy(s, p->host.gpu, p->dev, 48)) return false;
    p->sync    = md_gpu_stream_record(s);
    p->x       = x;
    p->y       = y;
    p->seq     = ++v->pick_seq;
    p->pending = true;
    return true;
}

bool md_render_view_pick_result(md_render_view_t v, md_render_pick_t* out) {
    if (!v || !out) return false;
    pick_slot_t* best = NULL;
    for (uint32_t i = 0; i < PICK_SLOTS; ++i) {
        pick_slot_t* p = &v->picks[i];
        if (!p->pending || !md_gpu_sync_is_complete(p->sync)) continue;
        p->pending = false;
        if (p->seq > v->pick_returned && (!best || p->seq > best->seq)) best = p;
    }
    if (!best) return false;
    v->pick_returned = best->seq;

    const uint32_t* w = (const uint32_t*)best->host.cpu;
    memset(out, 0, sizeof(*out));
    out->x    = best->x;
    out->y    = best->y;
    out->attr = UINT32_MAX;
    uint32_t slot1 = w[1] >> 20;
    if (slot1 == 0 || slot1 > best->meta_count) return true;   /* background */
    const pick_meta_t* m = &best->meta[slot1 - 1];
    uint64_t linear = ((uint64_t)(w[1] & 0x3FFFFu) << 32) | w[0];
    uint32_t cnt = m->count ? m->count : 1;
    out->key         = m->key;
    out->kind        = (md_render_kind_t)m->kind;
    out->batch_index = slot1 - 1;
    out->instance    = (uint32_t)(linear / cnt);
    out->element     = (uint32_t)(linear % cnt);
    out->id          = w[2];
    out->attr        = w[3];
    memcpy(&out->position[0], &w[4], sizeof(float));
    memcpy(&out->position[1], &w[5], sizeof(float));
    memcpy(&out->position[2], &w[6], sizeof(float));
    memcpy(&out->depth,       &w[7], sizeof(float));
    if (m->kind == MD_RENDER_EXTERNAL) {
        out->id   = out->element;
        out->attr = UINT32_MAX;
    }
    return true;
}

/* ---- Helpers ------------------------------------------------------------------------ */

static uint16_t f32_to_f16(float f) {
    uint32_t x;
    memcpy(&x, &f, 4);
    uint32_t sign = (x >> 16) & 0x8000u;
    int32_t  e    = (int32_t)((x >> 23) & 0xFF) - 127 + 15;
    uint32_t m    = x & 0x7FFFFFu;
    if (((x >> 23) & 0xFF) == 0xFF) return (uint16_t)(sign | 0x7C00u | (m ? 0x200u : 0u));
    if (e >= 31) return (uint16_t)(sign | 0x7C00u);
    if (e <= 0) {
        if (e < -10) return (uint16_t)sign;
        m |= 0x800000u;
        uint32_t shift = (uint32_t)(14 - e);
        uint32_t h = m >> shift;
        if ((m >> (shift - 1)) & 1u) h += 1;     /* round half up */
        return (uint16_t)(sign | h);
    }
    uint32_t h = sign | ((uint32_t)e << 10) | (m >> 13);
    if (m & 0x1000u) h += 1;                       /* round half up; carries into the exponent correctly */
    return (uint16_t)h;
}

void md_render_pack_velocity(const float v[3], uint32_t out[2]) {
    out[0] = (uint32_t)f32_to_f16(v[0]) | ((uint32_t)f32_to_f16(v[1]) << 16);
    out[1] = (uint32_t)f32_to_f16(v[2]);
}
