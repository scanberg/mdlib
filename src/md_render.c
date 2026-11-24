#include "md_render.h"
#include "md_render_debug.h"

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
    uint32_t      total_clusters, total_cluster_work, refit_mode;
    float         jitter_px;        /* how far this raster's jitter moves the image: max |jx|, |jy| */
    float         occl_margin;      /* phase 1: occluded only when nearer by this many radii too */
    float         occl_min_px;      /* elements no wider are not tested for occlusion (clusters are) */
    uint32_t      dbg_probe;      /* md_render_debug_raster_probe */
    uint32_t      spl_group;      /* md_render_debug_spline_group */
    uint32_t      spl_order;      /* md_render_debug_spline_order: buckets */
    uint32_t      _vpad6;
    uint32_t      _vpad3[2];
} gpu_view_t;
_Static_assert(sizeof(gpu_view_t) == 560, "ViewData mirror");

/* Occlusion settings, measured on the M4 at 10M atoms (camera moving):
   - Elements no wider than occl_min_px pixels are not tested one by one
     (their clusters are). Testing one costs about what drawing it does, and
     a tiny element rejected against the previous raster's pyramid leaves a
     hole in phase 1's: one uncovered pixel and the farthest-depth pyramid
     occludes nothing above it, so phase 2 let through everything behind.
     Far spacefill 5.4 -> 3.3 ms (phase 2's clusters 110k -> 1.5k); licorice
     14.8 -> 14.4 moving, 16.8 -> 17.6 still; near unchanged.
   - Phase 1's test passes what is nearer than the pyramid by fewer than
     occl_margin radii: elements just behind the front surface draw in phase
     1, for the same holes. Far 6.2 -> 5.4 ms, near 4.9 -> 4.5 (before the
     size threshold; beyond 2 radii nothing more).
   Settable for md_render_bench's experiments (declared there, not in the
   header). */
static float g_occl_margin = 2.0f;
static float g_occl_min_px = 2.0f;
void md_render_debug_occlusion(float margin_radii, float min_px) { g_occl_margin = margin_radii; g_occl_min_px = min_px; }

/* A probe of the rasters' costs (md_render_debug.h): 1, the spline and
   polygon fragments skip their hit and keep the proxy's depth (what the
   vertex stage and the fill cost alone). Wrong pictures. */
static uint32_t g_raster_probe = 0u;
void md_render_debug_raster_probe(uint32_t mode) { g_raster_probe = mode; }

/* Pieces per box of a spline part in the hardware raster: a fragment tests
   its box's pieces alone. 0: one box the part, whose fragments test every
   piece their ray may hit. By the device (spline_group_for): the bench's
   setting overrides it (not in the header; above 16, the device's again). */
static uint32_t g_spline_group = ~0u;
void md_render_debug_spline_group(uint32_t pieces) { g_spline_group = pieces; }

/* Nearest first (6f908220): on the 1060 1 piece took 27 to 37% off the
   others at 64 protein copies near (36 ms against 49 to 56; the 3080, in
   the cull's order: 21%); on the M4 4 pieces and the whole part tied, 1
   piece 4% behind; on the HD 530 4 pieces and the whole part tied, 1 piece
   4 to 12% behind at 64 copies mid, fit and the membrane. AMD and the rest:
   unmeasured, NVIDIA's. */
static uint32_t spline_group_for(uint32_t vendor_id) {
    if (g_spline_group <= 16u) return g_spline_group;
    return vendor_id == 0x106Bu || vendor_id == 0x8086u ? 4u : 1u;
}

/* Phase 1's spline boxes drawn nearest first, in this many buckets of view
   depth (spline_order, spline_bucket; the bench's setting, not in the
   header), each part's boxes by its box's centre. 0: the cull's order. On
   the M4 (c0955d44, 4 pieces a box) nearest first took 18 to 21% off the
   raster of 64 protein copies, 11 to 16% off the nucleosome's, 6 to 7% off
   the membrane's, about nothing off one copy's: the depth test turns away a
   fragment that writes depth against what the pass drew before it. 256
   buckets beat 64 by 2 to 5%, 64 beat 16 more; each box by its own centre
   or near face bought nothing; render passes of their own, nearest first,
   nothing either (both removed). On the HD 530 (6f908220, 4 pieces a box)
   it took 26 to 29% off the protein's 64 copies, 18 to 22% off the
   nucleosome's, 5 to 8% off the membrane's, up to 11% off one copy. On the
   1060 (1 piece a box) a wash: 10% off the 64-copy membrane near and fit,
   2% off the protein near, 4 to 9% on the protein's and the nucleosome's
   mid and fit. (On both the spline boxes drawn again over the final depth,
   in a pass after, cost as much as the first time, where the M4 shaded 39%
   of them again: that line does not tell what one pass turns away there.)
   On for all; the 3080 and AMD to measure. At most 256. */
static uint32_t g_spline_order = 256u;
void md_render_debug_spline_order(uint32_t buckets) { g_spline_order = buckets > 256u ? 256u : buckets; }
static uint32_t spline_order_on(void) { return g_spline_order != 0u; }

/* Occlusion culling by the view's own measure (occlusion_wanted; off, it
   always runs unless the view's debug options say otherwise). */
static bool g_occl_auto = true;
void md_render_debug_occlusion_auto(bool on) { g_occl_auto = on; }

enum {
    F_INDEX        = 1u << 0,
    F_VELOCITY     = 1u << 1,
    F_RADIUS       = 1u << 2,
    F_ATTR_INDEX   = 1u << 3,
    F_VISIBLE0     = 1u << 4,
    F_VISIBLE1     = 1u << 5,
    F_MARKS        = 1u << 6,
    F_INDEX_MAP    = 1u << 7,
    F_ID           = 1u << 8,
    F_RADIUS_SCALE = 1u << 9,
    F_COLOR        = 1u << 10,
    F_PALETTE_IDX  = 1u << 11,
    F_CLUSTERS     = 1u << 12,
    F_SHARED       = 1u << 13,
    F_RADIUS_TABLE = 1u << 14,
    F_FADE_OUT     = 1u << 15,
    F_FADE_IN      = 1u << 16,
    F_FADE_SHRINK  = 1u << 17,
    F_FADE_THIN    = 1u << 18,
    F_SPL_AXIS     = 1u << 19,
    F_SPL_PROFILE  = 1u << 20,
    F_SPL_THROUGH  = 1u << 21,
    F_DASHED       = 1u << 22,
    F_DISPLACED    = 1u << 23,
    F_SIDE         = 1u << 24,
    F_CIRCLE       = 1u << 25,
    F_FADE_KEY     = 1u << 26,
    F_ID_VISIBLE0  = 1u << 27,
    F_ID_VISIBLE1  = 1u << 28,
    F_SHAPED       = 1u << 29,
    F_FLAT         = 1u << 30,
    F_RADIUS_PX    = 1u << 31,
};

typedef struct gpu_batch_t {
    uint32_t      kind, count, instance_count, flags;
    uint32_t      stride, vel_stride, attr_count, color_mode;
    md_gpu_addr_t position, velocity, radius, attr_index, color, palette_index, palette;
    md_gpu_addr_t visible0, visible1, marks, index_map, index, id, radius_scale;
    md_gpu_addr_t transforms, transforms_prev;
    float         s0, s1;
    uint32_t      uniform_color, palette_count;
    md_gpu_float4 tint;
    float         sharpness, desaturate;
    uint32_t      work_offset, rad_stride;
    md_gpu_addr_t cluster_offset, bounds, cstate, radius_table;
    uint32_t      cluster_count, cluster_work_offset, cluster_base, radius_table_count;
    float         fade_out_a, fade_out_b, fade_in_a, fade_in_b;   /* the bands for this camera (fade_bands) */
    md_gpu_addr_t spl_axis, spl_profile;                          /* SPLINES */
    uint32_t      poly_n;                                         /* POLYGONS: corners */
    float         dash, gap, end_scale, end_add;                  /* CAPSULES dashed: the lengths, the spheres at the ends */
    float         cap_offset;                                     /* CAPSULES displaced: by this many radii */
    md_gpu_addr_t cap_side;                                       /* along this side (float4 per capsule), optional */
    float         poly_inset;                                     /* POLYGONS circles: the inset, */
    uint32_t      poly_dashes;                                    /* the dashes over each edge */
    float         cap_trim;                                       /* CAPSULES displaced: their ends pulled in */
    float         spl_scale;                                      /* SPLINES: of the profile */
    md_gpu_addr_t fade_key;                                       /* F_FADE_KEY: the point standing for each point */
    md_gpu_addr_t id_visible0, id_visible1;                       /* F_ID_VISIBLE0, 1: bitsets over the ids */
    float         cap_end, cap_px;                                /* CAPSULES shaped: the second radius, the radius in pixels */
} gpu_batch_t;
_Static_assert(sizeof(gpu_batch_t) == 368, "Batch mirror");

typedef struct gpu_args_t {
    md_gpu_addr_t        view;
    md_gpu_addr_t        batches;
    uint32_t             mode;      /* raster: 0 visibility, 1 mask. cull etc.: phase 1 or 2 */
    uint32_t             list;      /* raster: 0 hw0, 1 hw1, 2 st0, 3 st1, 4 hw2, 5 st2, 6 hw3, 7 st3 */
    md_gpu_sampled_tex_t id;
    md_gpu_storage_tex_t color, normal, motion, mask;
    uint32_t             width, height, px, py;
    md_gpu_addr_t        pick_out;
    md_gpu_addr_t        lists;
    md_gpu_sampled_tex_t src, hiz;
    md_gpu_storage_tex_t dst;
    uint32_t             src_mip, src_w, src_h, dst_w, dst_h, _pad0;
    md_gpu_storage_tex_t out_color, out_surface;
} gpu_args_t;
_Static_assert(sizeof(gpu_args_t) == 160, "Args mirror");

typedef struct gpu_lists_t {
    md_gpu_addr_t cvis, cocc, cst, hw0, hw1, st0, st1, small, occ, counters, vis64, hw2, st2, hw3, st3;
    md_gpu_addr_t hw2b;             /* per hw2 entry, its spline part's record: box, part, pieces' cuts: 20 float4 (SPLINE_RECORD) */
    md_gpu_addr_t hw2g;             /* boxes of groups of the parts' pieces: record << 4 | group, 16 at most a part */
    md_gpu_addr_t tiles;            /* the resolve's 8 x 8 tiles with spline pixels, then those with polygons */
    uint32_t      tiles_cap, _lpad; /* tiles of the target */
} gpu_lists_t;

/* counters[] (md_render.slang C_*) */
enum {
    C_HW1 = 0, C_HW2 = 2, C_ST = 4, C_SMALL1 = 6, C_SMALL2 = 7, C_OCC = 8,
    C_HIDDEN = 9, C_OUTSIDE = 10, C_OCCLUDED = 11,
    C_CVIS1 = 12, C_CVIS2 = 13, C_COCC = 14, C_CST = 15,
    C_CMD = 16, C_STCMD = 32, C_GRID_SMALL = 40, C_GRID_OCC = 48,
    C_GRID_CVIS = 52, C_GRID_COCC = 60, C_GRID_CST = 64, C_CL_VIS = 72, C_CL_OCC = 74,
    C_ICMD = 80,      /* md_gpu_draw_indexed_cmd_t: capsules, phase 1 and 2: 80..89 */
    C_STICMD = 90,    /* the state list's capsules: 90..94 */
    C_LOD = 95,       /* elements dropped by their batch's fade */
    C_HWS1 = 96, C_HWS2 = 97, C_STS = 98,   /* spline parts: hardware raster, phase 1 and 2; state */
    C_ICMDS = 100,    /* md_gpu_draw_indexed_cmd_t: spline parts, phase 1 and 2: 100..109 */
    C_STICMDS = 110,  /* the state list's spline parts: 110..114 */
    C_HWP1 = 116, C_HWP2 = 117, C_STP = 118,  /* polygons: hardware raster, phase 1 and 2; state */
    C_ICMDP = 120,    /* md_gpu_draw_indexed_cmd_t: polygons, phase 1 and 2: 120..129 */
    C_STICMDP = 130,  /* the state list's polygons: 130..134 */
    C_HWSG1 = 136, C_HWSG2 = 137,   /* boxes of the hardware spline parts' pieces (hw2g): phase 1 and 2 */
    C_GRID_ORD = 140,               /* spline_order's grid: 140..142 */
    C_ZRANGE = 316,                 /* phase 1's visible clusters' view depths, for the order: 316, 317 */
    C_BKT = 320, C_BKO = 576,       /* phase 1's spline boxes' depth buckets: counts and starts, 256 each */
    C_TILES = 832,                  /* the resolve's tiles with spline, with polygon pixels: 832, 833 */
    C_COUNT = 836,
};
#define IMPLICIT_CLUSTER 64u
#define CTR(v, i) ((v)->counters + (md_gpu_addr_t)(i) * 4u)

#define MAX_BATCHES 4095u
#define PICK_SLOTS  4u
#define STAT_SLOTS  3u
#define CULL_GROUP  64u
#define DEFAULT_SMALL_PIXELS 6.0f

/* A hook after each pass of a draw, for md_render_bench's per-pass profile
   (it syncs the stream in the hook). Not in the header: a debugging aid, set
   by md_render_debug_pass_hook, which the bench declares itself. */
static void (*g_pass_hook)(void* user, const char* pass);
static void* g_pass_user;

void md_render_debug_pass_hook(void (*hook)(void* user, const char* pass), void* user) {
    g_pass_hook = hook;
    g_pass_user = user;
}

#define PASS_DONE(name) do { if (g_pass_hook) g_pass_hook(g_pass_user, name); } while (0)

/* =========================================================================
   Context
   ========================================================================= */

struct md_render {
    md_gpu_device_t   dev;
    md_allocator_i*   alloc;
    md_gpu_pipeline_t sphere_vis, sphere_mask, capsule_vis, capsule_mask, spline_vis, spline_mask, polygon_vis, polygon_mask, merge;
    md_gpu_kernel_t   k_resolve, k_mask_visible, k_export, k_pick, k_cull, k_finalize, k_raster_small, k_hiz, k_refit, k_cull_clusters;
    /* the heavy twins (md_render.slang, kind_taken): run only for frames with splines or polygons */
    md_gpu_kernel_t   k_refit_heavy, k_cull_heavy, k_resolve_heavy, k_resolve_poly, k_spline_order;
    bool              atomic_int64;    /* the compute raster is possible */
    uint32_t          vendor_id;       /* the device's (spline_group_for) */
};

enum { PIPE_VIS, PIPE_MASK };
static md_gpu_pipeline_t make_pipeline(md_gpu_device_t dev, md_gpu_shader_t vs, md_gpu_shader_t fs, int kind, const char* label) {
    const bool mask = kind == PIPE_MASK;
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
        pd.color_count     = 2;                         /* the ids, and the heavy kinds' normals (VisOut) */
        pd.color[0].format = MD_GPU_FORMAT_RG32_UINT;
        pd.color[1].format = MD_GPU_FORMAT_RG16_FLOAT;
    }
    return md_gpu_pipeline_create(dev, &pd);
}

md_render_ctx_t md_render_create(md_gpu_device_t device, const md_render_desc_t* desc) {
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
    md_render_ctx_t r = md_alloc(alloc, sizeof(struct md_render));
    if (!r) return NULL;
    memset(r, 0, sizeof(*r));
    r->dev   = device;
    r->alloc = alloc;

    r->sphere_vis   = make_pipeline(device, md_shader_md_render_sphere_vs_shader(),  md_shader_md_render_sphere_fs_shader(),       PIPE_VIS,  "md_render.sphere");
    r->sphere_mask  = make_pipeline(device, md_shader_md_render_sphere_vs_shader(),  md_shader_md_render_sphere_mask_fs_shader(),  PIPE_MASK, "md_render.sphere_mask");
    r->capsule_vis  = make_pipeline(device, md_shader_md_render_capsule_vs_shader(), md_shader_md_render_capsule_fs_shader(),      PIPE_VIS,  "md_render.capsule");
    r->capsule_mask = make_pipeline(device, md_shader_md_render_capsule_vs_shader(), md_shader_md_render_capsule_mask_fs_shader(), PIPE_MASK, "md_render.capsule_mask");
    r->spline_vis   = make_pipeline(device, md_shader_md_render_spline_vs_shader(),  md_shader_md_render_spline_fs_shader(),       PIPE_VIS,  "md_render.spline");
    r->spline_mask  = make_pipeline(device, md_shader_md_render_spline_vs_shader(),  md_shader_md_render_spline_mask_fs_shader(),  PIPE_MASK, "md_render.spline_mask");
    r->polygon_vis  = make_pipeline(device, md_shader_md_render_polygon_vs_shader(), md_shader_md_render_polygon_fs_shader(),      PIPE_VIS,  "md_render.polygon");
    r->polygon_mask = make_pipeline(device, md_shader_md_render_polygon_vs_shader(), md_shader_md_render_polygon_mask_fs_shader(), PIPE_MASK, "md_render.polygon_mask");
    r->merge        = make_pipeline(device, md_shader_md_render_merge_vs_shader(),   md_shader_md_render_merge_fs_shader(),        PIPE_VIS,  "md_render.merge");
    r->atomic_int64 = info.atomic_int64_min_max;
    r->vendor_id    = info.vendor_id;

    md_gpu_kernel_desc_t kd;
    kd = md_shader_md_render_resolve_kernel();      r->k_resolve      = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_mask_visible_kernel(); r->k_mask_visible = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_export_targets_kernel(); r->k_export     = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_pick_kernel();         r->k_pick         = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_cull_kernel();         r->k_cull         = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_finalize_kernel();     r->k_finalize     = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_hiz_reduce_kernel();   r->k_hiz          = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_refit_kernel();        r->k_refit        = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_cull_clusters_kernel(); r->k_cull_clusters = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_refit_heavy_kernel();   r->k_refit_heavy   = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_cull_heavy_kernel();    r->k_cull_heavy    = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_resolve_heavy_kernel(); r->k_resolve_heavy = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_resolve_poly_kernel();  r->k_resolve_poly  = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_spline_order_kernel();  r->k_spline_order  = md_gpu_kernel_create(device, &kd);
    if (r->atomic_int64) {
        kd = md_shader_md_render_raster_small_kernel(); r->k_raster_small = md_gpu_kernel_create(device, &kd);
        if (!r->k_raster_small) r->atomic_int64 = false;
    }

    if (!r->sphere_vis || !r->sphere_mask || !r->capsule_vis || !r->capsule_mask || !r->spline_vis || !r->spline_mask ||
        !r->polygon_vis || !r->polygon_mask || !r->merge ||
        !r->k_resolve || !r->k_mask_visible || !r->k_export || !r->k_pick || !r->k_cull || !r->k_finalize || !r->k_hiz || !r->k_refit || !r->k_cull_clusters ||
        !r->k_refit_heavy || !r->k_cull_heavy || !r->k_resolve_heavy || !r->k_resolve_poly || !r->k_spline_order) {
        MD_LOG_ERROR("md_render: pipeline or kernel creation failed: %s", md_gpu_last_error());
        md_render_destroy(r);
        return NULL;
    }
    return r;
}

void md_render_destroy(md_render_ctx_t r) {
    if (!r) return;
    if (r->sphere_vis)     md_gpu_pipeline_destroy(r->sphere_vis);
    if (r->sphere_mask)    md_gpu_pipeline_destroy(r->sphere_mask);
    if (r->capsule_vis)    md_gpu_pipeline_destroy(r->capsule_vis);
    if (r->capsule_mask)   md_gpu_pipeline_destroy(r->capsule_mask);
    if (r->spline_vis)     md_gpu_pipeline_destroy(r->spline_vis);
    if (r->spline_mask)    md_gpu_pipeline_destroy(r->spline_mask);
    if (r->polygon_vis)    md_gpu_pipeline_destroy(r->polygon_vis);
    if (r->polygon_mask)   md_gpu_pipeline_destroy(r->polygon_mask);
    if (r->k_refit_heavy)    md_gpu_kernel_destroy(r->k_refit_heavy);
    if (r->k_cull_heavy)     md_gpu_kernel_destroy(r->k_cull_heavy);
    if (r->k_resolve_heavy)  md_gpu_kernel_destroy(r->k_resolve_heavy);
    if (r->k_resolve_poly)   md_gpu_kernel_destroy(r->k_resolve_poly);
    if (r->k_spline_order)   md_gpu_kernel_destroy(r->k_spline_order);
    if (r->k_resolve)      md_gpu_kernel_destroy(r->k_resolve);
    if (r->k_mask_visible) md_gpu_kernel_destroy(r->k_mask_visible);
    if (r->k_export)       md_gpu_kernel_destroy(r->k_export);
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
    md_render_ctx_t     r;
    char                label[64];

    uint32_t            width, height;
    md_render_targets_t t;

    /* md_render_debug_view's: MD_RENDER_DEBUG_*, the compute raster's extent;
       a change rasterises at the next draw */
    uint32_t            dbg_options;
    float               dbg_small_px;
    bool                dbg_changed;

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
    md_gpu_addr_t       zero;              /* 256 zero bytes: the absent arrays (no_null) */
    md_gpu_addr_t       box_idx;           /* uint16 [36]: a capsule's box over its 8 corners */
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
        bool            occl_tested;   /* occlusion was on, against a pyramid (occlusion_wanted) */
    } stats[STAT_SLOTS];
    uint32_t            stats_returned;
    uint32_t            occl_seen;         /* the last raster occlusion_wanted read the counts of */
    uint32_t            occl_off;          /* rasters left without occlusion */
    uint64_t            occl_drawn;        /* elements drawn when it was turned off */

    pick_slot_t         picks[PICK_SLOTS];
    uint64_t            pick_seq;
    uint64_t            pick_returned;
};

md_render_view_t md_render_view_create(md_render_ctx_t r, const char* label) {
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

void md_render_debug_view(md_render_view_t v, uint32_t options, float small_pixels) {
    if (!v) return;
    v->dbg_changed |= options != v->dbg_options || small_pixels != v->dbg_small_px;
    v->dbg_options  = options;
    v->dbg_small_px = small_pixels;
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
    if (v->zero)      md_gpu_free(s, v->zero);
    if (v->box_idx)   md_gpu_free(s, v->box_idx);
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
    v->t.normal = make_target(s, MD_GPU_FORMAT_RG16_FLOAT,  RT | SA | ST, w, h, "md_render.normal");
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
    if (b->kind == MD_RENDER_SPHERES && b->style.spheres.radius_scale != 0.0f && b->points.radius) {   /* as fill_batch reads them */
        out->points.radius             = b->points.radius;
        out->points.radius_stride      = b->points.radius_stride;
        out->points.radius_table       = b->points.radius_table;
        out->points.radius_table_count = b->points.radius_table ? b->points.radius_table_count : 0;
    }
    out->points.attr_count      = b->points.attr_count;
    out->points.attr_index      = b->points.attr_index;
    out->points.visible[0]      = b->points.visible[0];
    out->points.visible[1]      = b->points.visible[1];
    out->points.marks           = b->points.marks;      /* an address: the mask pass draws from it */
    out->index           = b->index;
    out->index_map       = b->index_map;
    out->id              = b->id;
    out->id_visible[0]   = b->id_visible[0];
    out->id_visible[1]   = b->id_visible[1];
    out->cluster_offset  = b->cluster_offset;
    out->cluster_count   = b->cluster_offset ? b->cluster_count : 0;
    if (b->kind == MD_RENDER_SPHERES) {
        out->style.spheres.radius_scale = b->style.spheres.radius_scale;
        out->style.spheres.radius_add   = b->style.spheres.radius_add;
    } else if (b->kind == MD_RENDER_CAPSULES) {
        const md_render_capsule_style_t* c = &b->style.capsules;
        out->style.capsules.radius       = c->radius;
        out->style.capsules.radius_scale = c->radius_scale;
        out->style.capsules.dash_length  = c->dash_length;
        out->style.capsules.gap_length   = c->gap_length;
        out->style.capsules.dash_ends    = c->dash_ends;
        out->style.capsules.offset       = c->offset;
        out->style.capsules.side         = c->offset != 0.0f ? c->side : 0;
        out->style.capsules.trim         = c->offset != 0.0f ? c->trim : 0.0f;
        out->style.capsules.flat         = c->flat;
        out->style.capsules.tapered      = c->tapered;
        out->style.capsules.end_radius   = c->tapered ? c->end_radius : 0.0f;
        out->style.capsules.radius_px    = c->radius_px;
    } else if (b->kind == MD_RENDER_SPLINES) {
        out->style.splines = b->style.splines;
    } else if (b->kind == MD_RENDER_POLYGONS) {
        out->style.polygons = b->style.polygons;
    }
    /* not the instances: the view keeps them apart from its copy of the
       batch (remember), and compare checks their count and transforms */
    if (b->fade.mode) out->fade = b->fade;
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

/* The kinds md_render culls and rasterises; and which hardware raster draws
   them: 0 spheres (quads), 1 capsules (boxes), 2 spline parts and 3
   polygons (boxes, their own stages, in the heavy twins: md_render.slang,
   kind_taken) */
#define HW_KINDS 4
static bool culled_kind(md_render_kind_t k) {
    return k == MD_RENDER_SPHERES || k == MD_RENDER_CAPSULES || k == MD_RENDER_POLYGONS || k == MD_RENDER_SPLINES;
}
static uint32_t hw_kind(md_render_kind_t k) {
    return k == MD_RENDER_SPHERES ? 0u : k == MD_RENDER_SPLINES ? 2u : k == MD_RENDER_POLYGONS ? 3u : 1u;
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
    if (f->jitter == MD_RENDER_JITTER_GIVEN && !(isfinite(f->jitter_offset[0]) && isfinite(f->jitter_offset[1]))) {
        MD_LOG_ERROR("md_render: a given jitter offset is not finite");
        return false;
    }
    uint64_t total = 0;
    for (uint32_t i = 0; i < f->batch_count; ++i) total += culled_kind(f->batches[i].kind)
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
        if (b->fade.mode) {
            const float* bands[2] = { b->fade.in, b->fade.out };
            for (int k = 0; k < 2; ++k) {
                const float* r = bands[k];
                if (r[0] == 0.0f && r[1] == 0.0f) continue;
                if (!(r[0] > r[1] && r[1] > 0.0f)) {
                    MD_LOG_ERROR("md_render: batch %u: a fade band is {start, end} with start > end > 0 (got {%g, %g})", i, r[0], r[1]);
                    return false;
                }
            }
            if (b->fade.in[0] > 0.0f && b->fade.out[0] > 0.0f && b->fade.in[1] < b->fade.out[0]) {
                MD_LOG_ERROR("md_render: batch %u: its fade's `in` band must end where `out` begins or nearer", i);
                return false;
            }
            if (b->fade.key && b->kind != MD_RENDER_SPHERES && b->kind != MD_RENDER_CAPSULES) {
                MD_LOG_ERROR("md_render: batch %u: a fade's key is for spheres and capsules", i);
                return false;
            }
        }
        switch (b->kind) {
        case MD_RENDER_SPHERES:
        case MD_RENDER_CAPSULES:
            break;
        case MD_RENDER_SPLINES:
            if (b->count && !b->index) {
                MD_LOG_ERROR("md_render: batch %u: splines need an index (control points and their flags)", i);
                return false;
            }
            if (b->index_map || b->id_visible[0] || b->id_visible[1]) {
                MD_LOG_ERROR("md_render: batch %u: splines take no index map and no visibility by id", i);
                return false;
            }
            break;
        case MD_RENDER_POLYGONS:
            if (b->style.polygons.corners < 3 || b->style.polygons.corners > 8) {
                MD_LOG_ERROR("md_render: batch %u: polygons have 3 to 8 corners (got %u)", i, b->style.polygons.corners);
                return false;
            }
            if (!(b->style.polygons.inset >= 0.0f) || !isfinite(b->style.polygons.inset) || b->style.polygons.dashes > 16) {
                MD_LOG_ERROR("md_render: batch %u: a circle's inset is finite and not negative, its dashes 16 at most", i);
                return false;
            }
            break;
        case MD_RENDER_EXTERNAL:
            continue;
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
            if (!isfinite(c->offset) || !(c->trim >= 0.0f) || !(c->trim <= 1.0f)) {
                MD_LOG_ERROR("md_render: batch %u: the offset is finite, the trim in [0, 1]", i);
                return false;
            }
            if (c->dash_length < 0.0f || c->gap_length < 0.0f) {
                MD_LOG_ERROR("md_render: batch %u: dash and gap lengths are not negative", i);
                return false;
            }
            const bool shaped = c->flat || c->tapered || c->radius_px != 0.0f;
            if (shaped && (c->dash_length > 0.0f || c->offset != 0.0f)) {
                MD_LOG_ERROR("md_render: batch %u: shaped capsules (flat, tapered, radius_px) are neither dashed nor displaced", i);
                return false;
            }
            if (!(c->radius_px >= 0.0f && isfinite(c->radius_px)) || (c->tapered && !(c->end_radius >= 0.0f && isfinite(c->end_radius)))) {
                MD_LOG_ERROR("md_render: batch %u: radius_px and end_radius are finite and >= 0", i);
                return false;
            }
            if (c->tapered && c->radius_px != 0.0f) {
                MD_LOG_ERROR("md_render: batch %u: a capsule's radius in pixels is not tapered", i);
                return false;
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
    if (!culled_kind(b->kind)) return 0;
    return (uint64_t)b->count * instance_count(b);
}

/* Whether the heavy twins run: a batch of spline parts or polygons to draw */
static bool frame_has_heavy(const md_render_frame_t* f) {
    for (uint32_t i = 0; i < f->batch_count; ++i) {
        const md_render_kind_t k = f->batches[i].kind;
        if ((k == MD_RENDER_SPLINES || k == MD_RENDER_POLYGONS) && batch_work(&f->batches[i])) return true;
    }
    return false;
}

/* Store the frame as the one last drawn. */
/* Clusters of a batch: given, or runs of 64 consecutive elements */
static uint32_t batch_clusters(const md_render_batch_t* b) {
    if (!culled_kind(b->kind)) return 0;
    if (b->cluster_offset && b->cluster_count) return b->cluster_count;
    return (b->count + IMPLICIT_CLUSTER - 1) / IMPLICIT_CLUSTER;
}

/* Two batches whose cluster bounds are the same refit: the same points,
   index, clusters, visibility, state and style; not the same key, ids (unless
   they are what is visible), velocities, colours or instances */
static bool bounds_shared(const md_render_batch_t* a, const md_render_batch_t* b) {
    md_render_batch_t x, y;
    canonical_geometry(&x, a);
    canonical_geometry(&y, b);
    x.key = y.key = 0;
    if (!x.id_visible[0] && !x.id_visible[1]) x.id = 0;
    if (!y.id_visible[0] && !y.id_visible[1]) y.id = 0;
    x.points.velocity = y.points.velocity = 0;
    x.points.velocity_stride = y.points.velocity_stride = 0;
    memset(&x.fade, 0, sizeof(x.fade));                     /* the bounds do not depend on it */
    memset(&y.fade, 0, sizeof(y.fade));
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
    if (f->jitter == MD_RENDER_JITTER_GIVEN &&
        (l->jitter != MD_RENDER_JITTER_GIVEN || memcmp(f->jitter_offset, l->jitter_offset, sizeof(f->jitter_offset)) != 0)) return CHANGE_RASTER;
    if (f->batch_count != l->batch_count) return CHANGE_RASTER;
    if (f->content[MD_RENDER_CONTENT_GEOMETRY] != l->content[MD_RENDER_CONTENT_GEOMETRY]) return CHANGE_RASTER;
    if (v->dbg_changed) return CHANGE_RASTER;
    change_t c = CHANGE_NONE;
    for (uint32_t i = 0; i < f->batch_count; ++i) {
        md_render_batch_t x, y;
        canonical_geometry(&x, &f->batches[i]);
        canonical_geometry(&y, &v->batches[i]);
        if (memcmp(&x, &y, sizeof(x)) != 0) return CHANGE_RASTER;
        uint32_t n = instance_count(&f->batches[i]);
        if (16u * n != v->xform_offset[i + 1] - v->xform_offset[i]) return CHANGE_RASTER;
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
        memset(&x.fade, 0, sizeof(x.fade));                     /* a raster, not a refit */
        memset(&y.fade, 0, sizeof(y.fade));
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

/* No array the shaders read is null: an absent one points at a block of zeros,
   and the shaders gate its index by its flag to 0 (md_render.slang, gated).
   A load in a branch not taken then reads that block. lavapipe executes such
   loads, at the address of the group's first lane whatever its mask: the pick
   read attr_index[34] where there was none, at 0x88, and crashed. A driver
   that speculates loads out of branches would do the same. */
static void no_null(gpu_batch_t* g, md_gpu_addr_t zero) {
    md_gpu_addr_t* opt[] = { &g->position, &g->velocity, &g->radius, &g->attr_index, &g->color, &g->palette_index, &g->palette,
                             &g->visible0, &g->visible1, &g->marks, &g->index_map, &g->index, &g->id, &g->id_visible0, &g->id_visible1, &g->radius_scale,
                             &g->radius_table, &g->cluster_offset, &g->spl_axis, &g->spl_profile, &g->cap_side, &g->fade_key };
    for (size_t i = 0; i < sizeof(opt) / sizeof(opt[0]); ++i) if (!*opt[i]) *opt[i] = zero;
}

static void fill_batch(gpu_batch_t* g, const md_render_batch_t* b, uint32_t n_inst, md_gpu_addr_t xf, md_gpu_addr_t xf_prev) {
    memset(g, 0, sizeof(*g));
    const md_render_points_t* p = &b->points;
    g->kind           = (uint32_t)b->kind;
    g->count          = b->count;
    g->instance_count = n_inst;
    g->stride         = p->stride ? p->stride : 12u;
    g->vel_stride     = p->velocity_stride ? p->velocity_stride : 8u;
    g->rad_stride     = p->radius_stride ? p->radius_stride : 4u;
    g->attr_count     = p->attr_count ? p->attr_count : p->count;
    g->color_mode     = (uint32_t)b->color.mode;
    g->position       = p->position;
    g->velocity       = p->velocity;
    /* Spheres read point radii, but not with a radius_scale of 0 (licorice:
       radius_add alone); so do dashed capsules, for the spheres at their
       ends, likewise */
    const bool radii  = p->radius && ((b->kind == MD_RENDER_SPHERES && b->style.spheres.radius_scale != 0.0f) ||
                                      (b->kind == MD_RENDER_CAPSULES && b->style.capsules.dash_length > 0.0f &&
                                       b->style.capsules.dash_ends.radius_scale != 0.0f));
    g->radius         = radii ? p->radius : 0;
    g->radius_table   = radii ? p->radius_table : 0;
    g->radius_table_count = g->radius_table ? (p->radius_table_count < 0x10000u ? p->radius_table_count : 0x10000u) : 0;
    g->attr_index     = p->attr_index;
    g->color          = p->color;
    g->palette_index  = p->palette_index;
    g->palette        = b->color.palette;
    g->visible0       = p->visible[0];
    g->visible1       = p->visible[1];
    g->marks          = p->marks;
    g->index          = b->index;
    g->index_map      = b->index_map;
    g->id             = b->id;
    g->id_visible0    = b->id_visible[0];
    g->id_visible1    = b->id_visible[1];
    g->transforms      = xf;
    g->transforms_prev = xf_prev;
    uint32_t f = 0;
    if (b->index)          f |= F_INDEX;
    if (b->index_map)      f |= F_INDEX_MAP;
    if (b->id_visible[0])  f |= F_ID_VISIBLE0;
    if (b->id_visible[1])  f |= F_ID_VISIBLE1;
    if (p->velocity)       f |= F_VELOCITY;
    if (g->radius)         f |= F_RADIUS;
    if (g->radius_table)   f |= F_RADIUS_TABLE;
    if (p->attr_index)     f |= F_ATTR_INDEX;
    if (p->visible[0])     f |= F_VISIBLE0;
    if (p->visible[1])     f |= F_VISIBLE1;
    if (p->marks)          f |= F_MARKS;
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
        if (b->style.capsules.dash_length > 0.0f) {
            f |= F_DASHED;
            g->dash      = b->style.capsules.dash_length;
            g->gap       = b->style.capsules.gap_length;
            g->end_scale = b->style.capsules.dash_ends.radius_scale;
            g->end_add   = b->style.capsules.dash_ends.radius_add;
        }
        if (b->style.capsules.offset != 0.0f) {
            f |= F_DISPLACED;
            g->cap_offset = b->style.capsules.offset;
            g->cap_side   = b->style.capsules.side;
            g->cap_trim   = b->style.capsules.trim;
            if (g->cap_side) f |= F_SIDE;
        }
        const md_render_capsule_style_t* cs = &b->style.capsules;
        if (cs->flat || cs->tapered || cs->radius_px > 0.0f) {
            f |= F_SHAPED;
            if (cs->flat) f |= F_FLAT;
            g->cap_end = cs->tapered ? cs->end_radius : cs->radius;
            if (cs->radius_px > 0.0f) {
                f |= F_RADIUS_PX;
                g->cap_px = cs->radius_px;
            }
        }
    } else if (b->kind == MD_RENDER_SPLINES) {
        g->s0        = b->style.splines.radius;
        g->s1        = b->style.splines.flat;
        g->spl_scale = b->style.splines.scale > 0.0f ? b->style.splines.scale : 1.0f;
        if (b->style.splines.through) f |= F_SPL_THROUGH;
        g->spl_axis  = b->style.splines.axis;
        g->spl_profile = b->style.splines.profile;
        if (g->spl_axis)    f |= F_SPL_AXIS;
        if (g->spl_profile) f |= F_SPL_PROFILE;
    } else if (b->kind == MD_RENDER_POLYGONS) {
        g->s0     = b->style.polygons.radius;
        g->s1     = b->style.polygons.plate > 0.0f ? MD_MIN(b->style.polygons.plate, 1.0f) : 1.0f;   /* 0: flat */
        g->poly_n = b->style.polygons.corners;
        if (b->style.polygons.inset > 0.0f) {
            f |= F_CIRCLE;
            g->poly_inset  = b->style.polygons.inset;
            g->poly_dashes = b->style.polygons.dashes;
        }
    }
    if (b->fade.mode) {
        if (b->fade.out[0] > 0.0f) f |= F_FADE_OUT;
        if (b->fade.in[0]  > 0.0f) f |= F_FADE_IN;
        if (f & (F_FADE_OUT | F_FADE_IN)) {
            if (b->fade.mode & MD_RENDER_FADE_SHRINK) f |= F_FADE_SHRINK;
            if (b->fade.mode & MD_RENDER_FADE_THIN)   f |= F_FADE_THIN;
            if (b->fade.key) { f |= F_FADE_KEY; g->fade_key = b->fade.key; }
        }
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

/* A batch's fade bands for this camera. The scale (pixels per unit) at
   view-space depth d is px / d, px = |proj[5]| * height / 2 (unjittered), so
   a band {s0, s1} is the depths px / s0 .. px / s1, and how far through it
   an element at depth d is, t = saturate(a d + b), is linear in d.
   Orthographic: the scale is px everywhere, every element equally far
   through: a = 0, b = t. Both sides of a handover get the same a and b from
   the same two numbers. */
static void fade_band(float out_ab[2], const float band[2], float px, bool ortho) {
    const double d0 = 1.0 / band[0], d1 = 1.0 / band[1];     /* depths over px */
    if (ortho) {
        out_ab[0] = 0.0f;
        out_ab[1] = (float)((1.0 / px - d0) / (d1 - d0));
        return;
    }
    const double a = 1.0 / (px * (d1 - d0));
    out_ab[0] = (float)a;
    out_ab[1] = (float)(-px * d0 * a);
}

static void fade_bands(gpu_batch_t* g, const md_render_batch_t* b, const md_render_frame_t* f) {
    g->fade_out_a = g->fade_out_b = g->fade_in_a = g->fade_in_b = 0.0f;
    if (!(g->flags & (F_FADE_OUT | F_FADE_IN))) return;
    const float* P  = f->camera.proj;
    const float  px = fabsf(P[5]) * 0.5f * (float)f->height;
    const bool ortho = P[3] == 0.0f && P[7] == 0.0f && P[11] == 0.0f && P[15] == 1.0f;
    float ab[2];
    if (g->flags & F_FADE_OUT) { fade_band(ab, b->fade.out, px, ortho); g->fade_out_a = ab[0]; g->fade_out_b = ab[1]; }
    if (g->flags & F_FADE_IN)  { fade_band(ab, b->fade.in,  px, ortho); g->fade_in_a  = ab[0]; g->fade_in_b  = ab[1]; }
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
    if (!v->zero) {
        v->zero = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, 256).gpu;
        if (!v->zero || !md_gpu_memset(s, v->zero, 0, 256)) return false;
    }
    uint32_t work = 0, cbase = 0, cwork = 0;
    for (uint32_t i = 0; i < f->batch_count; ++i) {
        const md_render_batch_t* b = &f->batches[i];
        md_gpu_addr_t xf = v->table[i].transforms, xp = v->table[i].transforms_prev;
        fill_batch(&v->table[i], b, instance_count(b), xf, xp);
        fade_bands(&v->table[i], b, f);
        gpu_batch_t* g = &v->table[i];
        g->work_offset         = work;
        g->cluster_count       = batch_clusters(b);
        g->cluster_offset      = (b->cluster_offset && b->cluster_count) ? b->cluster_offset : 0;
        if (g->cluster_offset) g->flags |= F_CLUSTERS;
        no_null(g, v->zero);
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
   its kind (hardware, state, small, occluded), 4 bytes an entry: 16 bytes per
   element in all (8 for polygons: no compute raster, no retest; 392 for
   spline parts, whose hardware entries carry their box, the part itself and
   its pieces' cuts of the box, 320 bytes, and up to 16 boxes of their
   pieces, 4 bytes each), plus 24 per cluster. */
static bool ensure_lists(md_gpu_stream_t s, md_render_view_t v, const uint64_t work[HW_KINDS], uint64_t cwork, bool compute) {
    const uint64_t tot = work[0] + work[1];
    /* Cluster entries are 8 bytes, element entries 4 (a global work index) */
    const uint64_t tiles = (uint64_t)((v->width + 7u) / 8u) * ((v->height + 7u) / 8u);
    const uint64_t sizes[16] = { cwork, cwork, cwork, work[0], work[1], work[0], work[1], tot, tot, work[2], work[2], work[3], work[3], work[2], 16u * work[2],
                                 (work[2] || work[3]) ? 2u * tiles : 0u };
    const uint32_t esize[16] = { 8, 8, 8, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 320, 4, 4 };   /* SPLINE_RECORD float4; spline_pieces' 16 */
    size_t need = 0;
    for (int i = 0; i < 16; ++i) need += (size_t)align_addr((sizes[i] + 1) * esize[i]);
    if (!gpu_realloc(s, &v->list_mem, &v->list_mem_cap, need)) return false;
    gpu_lists_t l;
    md_gpu_addr_t o = v->list_mem;
    md_gpu_addr_t* dst[16] = { &l.cvis, &l.cocc, &l.cst, &l.hw0, &l.hw1, &l.st0, &l.st1, &l.small, &l.occ, &l.hw2, &l.st2, &l.hw3, &l.st3, &l.hw2b, &l.hw2g, &l.tiles };
    for (int i = 0; i < 16; ++i) { *dst[i] = o; o += align_addr((sizes[i] + 1) * esize[i]); }
    l.tiles_cap = (uint32_t)sizes[15];
    l._lpad     = 0;
    if (!v->counters) v->counters = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, C_COUNT * 4u).gpu;
    if (!v->lists_buf) v->lists_buf = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, sizeof(gpu_lists_t)).gpu;
    if (!v->counters || !v->lists_buf) return false;
    if (!v->box_idx) {
        /* 12 triangles counter-clockwise seen from outside; corner c is
           (c & 1, c >> 1 & 1, c >> 2 & 1) in the capsule's frame (u, v, axis).
           Indexed, the vertex stage runs once per corner, not per index. */
        static const uint16_t box[36] = {
            0, 4, 6,  0, 6, 2,   1, 3, 7,  1, 7, 5,   0, 1, 5,  0, 5, 4,
            2, 6, 7,  2, 7, 3,   0, 2, 3,  0, 3, 1,   4, 5, 7,  4, 7, 6,
        };
        v->box_idx = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, sizeof(box)).gpu;
        if (!v->box_idx || !md_gpu_upload(s, v->box_idx, box, sizeof(box))) return false;
    }
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

/* The refit: eight clusters a group, refit_groups(n) of the shader */
static bool pass_per_cluster(md_gpu_stream_t s, md_render_view_t v, md_gpu_kernel_t k, uint32_t mode, uint32_t clusters) {
    gpu_args_t a = base_args(v);
    a.mode = mode;
    uint32_t groups = (clusters + 7u) / 8u;
    if (groups > 65535u) groups = 65535u;
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
static bool run_state_lists(md_gpu_stream_t s, md_render_view_t v, uint32_t cluster_work, bool splines) {
    md_render_ctx_t r = v->r;
    return pass(s, v, r->k_cull_clusters, 3, cluster_work, 0) && fin(s, v, 31) &&
           pass(s, v, r->k_cull, 3, 0, C_GRID_CST) && (!splines || pass(s, v, r->k_cull_heavy, 3, 0, C_GRID_CST)) && fin(s, v, 32);
}

static bool run_raster_small(md_gpu_stream_t s, md_render_view_t v, uint32_t phase) {
    gpu_args_t a = base_args(v);
    a.mode = phase;
    return md_gpu_launch_indirect(s, v->r->k_raster_small, CTR(v, C_GRID_SMALL + (phase - 1) * 4), &a, sizeof(a));
}

/* Phase 1 clears; phase 2 adds. The compute raster's result goes in first, as
   depth and ids, then each kind in one indirect draw whose instance count the
   cull wrote. The CPU's cost is independent of the number of batches. */
static bool run_raster(md_gpu_stream_t s, md_render_view_t v, uint32_t phase, bool compute, const bool has_kind[HW_KINDS]) {
    md_render_ctx_t r = v->r;
    md_gpu_render_desc_t rd = {0};
    rd.color_count       = 2;
    rd.color[0].texture  = v->t.id;
    rd.color[0].load     = phase == 1 ? MD_GPU_LOAD_CLEAR : MD_GPU_LOAD;
    rd.color[0].store    = MD_GPU_STORE;
    rd.color[1].texture  = v->t.normal;                 /* the heavy kinds' hit normals (VisOut); the resolve's for the rest */
    rd.color[1].load     = phase == 1 ? MD_GPU_LOAD_CLEAR : MD_GPU_LOAD;
    rd.color[1].store    = MD_GPU_STORE;
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
    for (uint32_t k = 0; k < HW_KINDS; ++k) {
        if (!has_kind[k]) continue;
        ds.cull = k == 0 ? MD_GPU_CULL_NONE : MD_GPU_CULL_BACK;
        md_gpu_set_draw_state(s, &ds);
        a.mode = 0;
        a.list = k == 2 ? 4 : k == 3 ? 6 : k;
        if (k == 0)      md_gpu_draw_indirect(s, r->sphere_vis, CTR(v, C_CMD + (phase - 1) * 2 * 4), 1, &a, sizeof(a));
        else if (k == 1) md_gpu_draw_indexed_indirect(s, r->capsule_vis, v->box_idx, MD_GPU_INDEX_U16, CTR(v, C_ICMD + (phase - 1) * 5), 1, &a, sizeof(a));
        else if (k == 2) md_gpu_draw_indexed_indirect(s, r->spline_vis, v->box_idx, MD_GPU_INDEX_U16, CTR(v, C_ICMDS + (phase - 1) * 5), 1, &a, sizeof(a));
        else             md_gpu_draw_indexed_indirect(s, r->polygon_vis, v->box_idx, MD_GPU_INDEX_U16, CTR(v, C_ICMDP + (phase - 1) * 5), 1, &a, sizeof(a));
    }
    return md_gpu_render_end(s);
}

/* The main resolve, then its heavy twin over the spline and polygon pixels */
/* The light kinds' kernel, then a kernel for each heavy kind the frame has:
   splines and polygons together held 255 registers on the 1060 and spilled.
   The last pass is left for the caller to mark (md_render_debug_pass_hook). */
static const char* run_resolve(md_gpu_stream_t s, md_render_view_t v, const md_render_frame_t* f) {
    bool spl = false, poly = false;
    /* (the heavy kernels: the raster's depth, the point of the hit it took) */
    for (uint32_t i = 0; i < f->batch_count; ++i) {
        if (!batch_work(&f->batches[i])) continue;
        spl  |= f->batches[i].kind == MD_RENDER_SPLINES;
        poly |= f->batches[i].kind == MD_RENDER_POLYGONS;
    }
    gpu_args_t a = base_args(v);
    a.src = md_gpu_texture_sampled(v->t.depth);
    md_gpu_grid_t g = md_gpu_grid_for(v->r->k_resolve, v->width, v->height, 1);
    if (!md_gpu_launch(s, v->r->k_resolve, g, &a, sizeof(a))) return NULL;
    if (!spl && !poly) return "resolve";
    PASS_DONE("resolve");
    /* the heavy ones over the tiles the main one listed (C_TILES): a group a
       tile, so many groups taking them in turn (Args.mode) */
    const uint32_t tiles  = ((v->width + 7u) / 8u) * ((v->height + 7u) / 8u);
    const uint32_t groups = tiles < 4096u ? tiles : 4096u;
    a.mode = groups;
    if (spl) {
        if (!md_gpu_launch(s, v->r->k_resolve_heavy, md_gpu_grid(groups, 1, 1), &a, sizeof(a))) return NULL;
        if (!poly) return "resolve splines";
        PASS_DONE("resolve splines");
    }
    if (!md_gpu_launch(s, v->r->k_resolve_poly, md_gpu_grid(groups, 1, 1), &a, sizeof(a))) return NULL;
    return "resolve polygons";
}

/* The visible part from the ids; the hidden part by rasterising the state
   lists of the cull without a depth test, so it costs as much as the
   selection and not the system. */
static bool run_masks(md_gpu_stream_t s, md_render_view_t v, const md_render_frame_t* f) {
    gpu_args_t a = base_args(v);
    md_gpu_grid_t g = md_gpu_grid_for(v->r->k_mask_visible, v->width, v->height, 1);
    if (!md_gpu_launch(s, v->r->k_mask_visible, g, &a, sizeof(a))) return false;

    bool has[HW_KINDS] = { false, false, false, false };
    for (uint32_t i = 0; i < f->batch_count; ++i) {
        const md_render_batch_t* b = &f->batches[i];
        if (b->count && b->points.marks) {
            if (culled_kind(b->kind)) has[hw_kind(b->kind)] = true;
        }
    }
    if (!has[0] && !has[1] && !has[2] && !has[3]) return true;
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
    for (uint32_t k = 0; k < HW_KINDS; ++k) {
        if (!has[k]) continue;
        ds.cull = k == 0 ? MD_GPU_CULL_NONE : MD_GPU_CULL_BACK;
        md_gpu_set_draw_state(s, &ds);
        a.mode = 1;
        a.list = k == 2 ? 5 : k == 3 ? 7 : 2 + k;
        if (k == 0)      md_gpu_draw_indirect(s, v->r->sphere_mask, CTR(v, C_STCMD), 1, &a, sizeof(a));
        else if (k == 1) md_gpu_draw_indexed_indirect(s, v->r->capsule_mask, v->box_idx, MD_GPU_INDEX_U16, CTR(v, C_STICMD), 1, &a, sizeof(a));
        else if (k == 2) md_gpu_draw_indexed_indirect(s, v->r->spline_mask, v->box_idx, MD_GPU_INDEX_U16, CTR(v, C_STICMDS), 1, &a, sizeof(a));
        else             md_gpu_draw_indexed_indirect(s, v->r->polygon_mask, v->box_idx, MD_GPU_INDEX_U16, CTR(v, C_STICMDP), 1, &a, sizeof(a));
    }
    return md_gpu_render_end(s);
}

static void record_stats(md_gpu_stream_t s, md_render_view_t v, uint64_t elements, uint32_t clusters, bool refitted, bool occl_tested) {
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
    v->stats[best].occl_tested = occl_tested;
}

/* Whether this raster culls by occlusion. Its fixed cost -- two pyramids,
   phase 2's passes, its render pass (on the M4 the second raster pass alone
   loads and stores the targets, 0.3 ms at 1080p) -- is 0.3 to 0.7 ms on the
   M4 and the 1060 at 1080p, and its tests grow with the elements; what it
   culls can be worth less. Occlusion always against none at 6f908220, ms
   (M4 / 1060; the share of the elements it culled; their count weighed as
   below):
     cartoon, one copy and 64 (0 to 1%)                     +0.3..+0.5 / +0.3..+0.4
     the 64-copy membrane (16%; 13855)                      +0.5 / 0
     100k atoms, licorice (4 to 22%; 3k to 4k)              +0.5 / +0.5..+0.6
     100k, spacefill far (76%; 38k)                         +0.35 / +0.15
     100k, spacefill near (84%; 72k)                        0 / -1.7
     1M, licorice near (8%; 35k)                            +0.55 / +0.6
     1M, spacefill far (85%; 106k) / near (93%; 802k)       -0.15 / -1.0, -3.4 / -8.8
     10M, licorice near (12 to 17%; 180k to 290k)           0 / +1.1..+1.5
     10M, licorice into spacefill near (46 to 49%; 1M)      -4.9 / -8.2
     10M, spacefill (78 to 97%; 1M to 8.3M)                 -4.3..-30 / -11..-66
   So it pays where it culls a good share of the elements and many of them:
   at least OCCL_SHARE of those it culls or draws, and OCCL_LITTLE_PER_PX a pixel of
   the target (about 52k at 1080p: the fixed cost grows with the pixels),
   the culled weighed as the drawn ones go, those of the compute raster at
   an eighth (their raster is that much cheaper). A raster that tested it
   and found less turns it off for the next OCCL_OFF rasters; then it is
   tested again, sooner when the elements drawn more than double. (The
   first takes, 10k culled a frame and then 0.05 a pixel without the share,
   kept it on for the 10M licorice that the 1060 drew 1.5 ms faster without,
   and off for the 100k spacefill near that it drew 1.7 ms faster with.)
   The image does not depend on it. The counts are the stats' (record_stats),
   read once complete, never waited for. */
#define OCCL_SHARE         0.3
#define OCCL_LITTLE_PER_PX 0.025
#define OCCL_OFF           60u
static bool occlusion_wanted(md_render_view_t v) {
    if (!g_occl_auto) return true;
    int best = -1;
    for (uint32_t i = 0; i < STAT_SLOTS; ++i) {
        if (!v->stats[i].host.cpu || v->stats[i].raster <= v->occl_seen || !md_gpu_sync_is_complete(v->stats[i].sync)) continue;
        if (best < 0 || v->stats[i].raster > v->stats[best].raster) best = (int)i;
    }
    if (best >= 0) {
        const uint32_t* c = (const uint32_t*)v->stats[best].host.cpu;
        const uint64_t hw    = (uint64_t)c[C_HW1] + c[C_HW1 + 1] + c[C_HWS1] + c[C_HWP1] + c[C_HW2] + c[C_HW2 + 1] + c[C_HWS2] + c[C_HWP2];
        const uint64_t small = (uint64_t)c[C_SMALL1] + c[C_SMALL2];
        const uint64_t drawn = hw + small;
        v->occl_seen = v->stats[best].raster;
        if (v->stats[best].occl_tested) {
            const double culled = (double)c[C_OCCLUDED];
            const double weight = drawn ? ((double)hw + (double)small / 8.0) / (double)drawn : 1.0;
            if (culled < OCCL_SHARE * (culled + (double)drawn) || culled * weight < OCCL_LITTLE_PER_PX * (double)v->width * (double)v->height) {
                v->occl_off   = OCCL_OFF;
                v->occl_drawn = drawn;
            }
        } else if (v->occl_off && (double)drawn > 2.0 * (double)v->occl_drawn + OCCL_LITTLE_PER_PX * (double)v->width * (double)v->height) {
            v->occl_off = 0;
        }
    }
    if (v->occl_off) { v->occl_off--; return false; }
    return true;
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
    out->faded       = c[C_LOD];
    out->occluded    = c[C_OCCLUDED];
    out->hardware[0] = c[C_HW1] + c[C_HW1 + 1] + c[C_HWS1] + c[C_HWP1];
    out->hardware[1] = c[C_HW2] + c[C_HW2 + 1] + c[C_HWS2] + c[C_HWP2];
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

    const bool occlusion = !(v->dbg_options & MD_RENDER_DEBUG_NO_OCCLUSION) && (!raster || occlusion_wanted(v));
    const bool compute   = !(v->dbg_options & MD_RENDER_DEBUG_NO_COMPUTE_RASTER) && v->r->atomic_int64;
    bool hiz_valid = false;
    if (resized) v->hiz_ok = false;

    if (raster) {
        /* Jitter, in pixels. Sample 0 is pixel centres, so a view drawn once
           (or a test) sees exactly what its camera says. */
        float jx = 0.0f, jy = 0.0f;
        if (f->jitter == MD_RENDER_JITTER_GIVEN) {
            jx = f->jitter_offset[0];
            jy = f->jitter_offset[1];
        } else if (f->jitter != MD_RENDER_JITTER_NONE && sample > 0) {
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
        gv.small_px        = v->dbg_small_px > 0.0f ? v->dbg_small_px : DEFAULT_SMALL_PIXELS;
        gv.jitter_px       = fabsf(jx) > fabsf(jy) ? fabsf(jx) : fabsf(jy);
        gv.occl_margin     = g_occl_margin;
        gv.occl_min_px     = g_occl_min_px;
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
        gv.dbg_probe  = g_raster_probe;
        gv.spl_group  = spline_group_for(v->r->vendor_id);
        gv.spl_order  = g_spline_order;
        md_gpu_upload(s, v->view_buf, &gv, sizeof(gv));

        uint64_t work[HW_KINDS] = { 0, 0, 0, 0 };
        bool has_kind[HW_KINDS] = { false, false, false, false };
        for (uint32_t i = 0; i < f->batch_count; ++i) {
            const md_render_batch_t* b = &f->batches[i];
            uint32_t k = hw_kind(b->kind);
            work[k] += batch_work(b);
            if (batch_work(b)) has_kind[k] = true;
        }
        const bool spl = has_kind[2] || has_kind[3];   /* the heavy twins run */
        const uint32_t total = (uint32_t)(work[0] + work[1] + work[2] + work[3]);
        const md_gpu_addr_t list_mem_before = v->list_mem;
        if (!ensure_lists(s, v, work, gv.total_cluster_work, compute)) return false;
        /* State lists: rebuilt only when the state bits or the geometry
           changed (they do not depend on the camera); their counts and
           commands survive the counters' reset otherwise */
        const bool state_lists = refit || state_changed || !v->state_lists_ok || v->list_mem != list_mem_before;
        if (state_lists) {
            md_gpu_memset(s, v->counters, 0, C_COUNT * 4u);
        } else {
            static const uint32_t keep[][2] = { { C_ST, C_ST + 2 }, { C_CST, C_CST + 1 }, { C_STCMD, C_STCMD + 8 }, { C_GRID_CST, C_GRID_CST + 3 },
                                                { C_STICMD, C_STICMD + 5 }, { C_STS, C_STS + 1 }, { C_STICMDS, C_STICMDS + 5 },
                                                { C_STP, C_STP + 1 }, { C_STICMDP, C_STICMDP + 5 } };
            uint32_t w = 0;
            for (size_t i = 0; i < sizeof(keep) / sizeof(keep[0]); ++i) {
                if (keep[i][0] > w) md_gpu_memset(s, CTR(v, w), 0, (keep[i][0] - w) * 4u);
                w = keep[i][1];
            }
            if (w < C_COUNT) md_gpu_memset(s, CTR(v, w), 0, (C_COUNT - w) * 4u);
        }
        if (compute) md_gpu_memset(s, v->vis64, 0, (size_t)v->width * v->height * 8u);
        PASS_DONE("setup");

        /* Cluster bounds: refitted when the geometry changed, else only their
           state bits when those did. Not when only the camera moved. */
        md_render_ctx_t r = v->r;
        if ((refit || state_changed) && !pass_per_cluster(s, v, r->k_refit, 0, gv.total_clusters)) return false;
        if (refit || state_changed) PASS_DONE("refit");
        if ((refit || state_changed) && spl && !pass_per_cluster(s, v, r->k_refit_heavy, 0, gv.total_clusters)) return false;
        if ((refit || state_changed) && spl) PASS_DONE("refit heavy");
        v->bounds_valid = true;

        /* Phase 1: clusters, then the elements of the visible ones, against
           the previous raster's pyramid (when there is one) */
        if (!pass(s, v, r->k_cull_clusters, 1, gv.total_cluster_work, 0) || !fin(s, v, 11)) return false;
        PASS_DONE("clusters1");
        if (!pass(s, v, r->k_cull, 1, 0, C_GRID_CVIS)) return false;
        if (spl) PASS_DONE("cull1");
        if (spl && !pass(s, v, r->k_cull_heavy, 1, 0, C_GRID_CVIS)) return false;
        if (!fin(s, v, 12)) return false;
        if (has_kind[2] && spline_order_on() && !pass(s, v, r->k_spline_order, 1, 0, C_GRID_ORD)) return false;
        PASS_DONE(spl ? "cull1 heavy" : "cull1");
        if (compute && !run_raster_small(s, v, 1)) return false;
        if (compute) PASS_DONE("small1");
        if (!run_raster(s, v, 1, compute, has_kind)) return false;
        PASS_DONE("raster1");

        if (occlusion) {
            if (!build_hiz(s, v)) return false;
            PASS_DONE("hiz1");
            if (hiz_valid) {
                /* Phase 2: the clusters and the elements phase 1 found
                   occluded, against this raster */
                if (!pass(s, v, r->k_cull_clusters, 2, 0, C_GRID_COCC) || !fin(s, v, 21)) return false;
                PASS_DONE("clusters2");
                if (!pass(s, v, r->k_cull, 2, 0, C_GRID_OCC)) return false;
                PASS_DONE("retest2");
                if (!pass(s, v, r->k_cull, 4, 0, C_GRID_CVIS + 4)) return false;
                if (spl) PASS_DONE("cull2");
                if (spl && !pass(s, v, r->k_cull_heavy, 4, 0, C_GRID_CVIS + 4)) return false;
                if (!fin(s, v, 22)) return false;
                PASS_DONE(spl ? "cull2 heavy" : "cull2");
                if (compute && !run_raster_small(s, v, 2)) return false;
                if (compute) PASS_DONE("small2");
                if (!run_raster(s, v, 2, compute, has_kind)) return false;
                PASS_DONE("raster2");
                if (!build_hiz(s, v)) return false;     /* for the next raster */
                PASS_DONE("hiz2");
            }
            v->hiz_ok = true;
        } else {
            v->hiz_ok = false;
        }
        if (state_lists) {
            if (!run_state_lists(s, v, gv.total_cluster_work, spl)) return false;
            v->state_lists_ok = true;
            PASS_DONE("state lists");
        }
        v->raster_count++;
        record_stats(s, v, total, gv.total_cluster_work, refit, hiz_valid);

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
        if (!pass_per_cluster(s, v, v->r->k_refit, 0, clusters)) return false;
        PASS_DONE("state refit");
        if (frame_has_heavy(f) && !pass_per_cluster(s, v, v->r->k_refit_heavy, 0, clusters)) return false;
        if (frame_has_heavy(f)) PASS_DONE("state refit heavy");
        md_gpu_memset(s, CTR(v, C_ST), 0, 8u);
        md_gpu_memset(s, CTR(v, C_STS), 0, 4u);
        md_gpu_memset(s, CTR(v, C_STP), 0, 4u);
        md_gpu_memset(s, CTR(v, C_CST), 0, 4u);
        if (!run_state_lists(s, v, cwork, frame_has_heavy(f))) return false;
        PASS_DONE("state lists");
    }

    if (resolve) {
        const char* last = run_resolve(s, v, f);
        if (!last) return false;
        PASS_DONE(last);
    }
    if (mask && !run_masks(s, v, f)) return false;
    if (mask) PASS_DONE("mask");

    if (!remember(v, f)) return false;
    v->has_raster = true;
    v->dbg_changed = false;
    v->sample     = sample;

    if (out) {
        out->work = (raster ? MD_RENDER_WORK_RASTER : 0u) | (resolve ? MD_RENDER_WORK_RESOLVE : 0u) | (mask ? MD_RENDER_WORK_MASK : 0u);
        out->sample    = sample;
        out->jitter[0] = v->jitter[0];
        out->jitter[1] = v->jitter[1];
    }
    return true;
}

/* ---- Marks ------------------------------------------------------------- */

void md_render_marks_write(uint8_t* dst, size_t count, md_bitset_t state0, md_bitset_t state1) {
    if (!dst) return;
    for (size_t base = 0; base < count; base += 64) {
        const uint64_t w0 = md_bitset_word(state0, (uint32_t)base);
        const uint64_t w1 = md_bitset_word(state1, (uint32_t)base);
        const size_t   n  = count - base < 64 ? count - base : 64;
        for (size_t i = 0; i < n; ++i) dst[base + i] = (uint8_t)(((w0 >> i) & 1u) | (((w1 >> i) & 1u) << 1));
    }
}

void md_render_marks_tag(uint8_t* dst, size_t count, md_bitset_t set, uint32_t tag) {
    if (!dst) return;
    const uint8_t t = MD_RENDER_MARK_TAG(tag), keep = (uint8_t)~MD_RENDER_MARK_TAG(7u);
    for (size_t base = 0; base < count; base += 64) {
        uint64_t w = md_bitset_word(set, (uint32_t)base);
        if (count - base < 64) w &= (1ull << (count - base)) - 1u;
        while (w) {
            const size_t i = base + md_ctz64(w);
            dst[i] = (uint8_t)((dst[i] & keep) | t);
            w &= w - 1u;
        }
    }
}

md_render_targets_t md_render_view_targets(md_render_view_t v) {
    md_render_targets_t t = {0};
    if (v) t = v->t;
    return t;
}

static bool export_dst_ok(const md_render_view_t v, md_gpu_texture_t t, md_gpu_format_t format, const char* what) {
    if (!t) return true;
    const md_gpu_texture_desc_t* d = md_gpu_texture_desc(t);
    if (!d || d->type != MD_GPU_TEX_2D || d->format != format || !(d->usage & MD_GPU_TEX_STORAGE) || d->width != v->width || d->height != v->height) {
        MD_LOG_ERROR("md_render_view_export: %s is to be a 2D STORAGE texture of %u x %u (R32_FLOAT, RGBA8_UNORM, RGBA16_FLOAT: depth, color, surface)",
                     what, v->width, v->height);
        return false;
    }
    return true;
}

bool md_render_view_export(md_gpu_stream_t s, md_render_view_t v, const md_render_export_t* dst) {
    if (!s || !v || !dst) return false;
    if (!v->has_raster) {
        MD_LOG_ERROR("md_render_view_export: the view has not been drawn");
        return false;
    }
    if (!export_dst_ok(v, dst->depth,   MD_GPU_FORMAT_R32_FLOAT,    "depth")  ||
        !export_dst_ok(v, dst->color,   MD_GPU_FORMAT_RGBA8_UNORM,  "color")  ||
        !export_dst_ok(v, dst->surface, MD_GPU_FORMAT_RGBA16_FLOAT, "surface")) return false;
    if (!dst->depth && !dst->color && !dst->surface) return true;
    gpu_args_t a = base_args(v);
    a.src  = md_gpu_texture_sampled(v->t.depth);
    a.mode = (dst->depth ? 1u : 0u) | (dst->color ? 2u : 0u) | (dst->surface ? 4u : 0u);
    if (dst->depth)   a.dst         = md_gpu_texture_storage(dst->depth, 0);
    if (dst->color)   a.out_color   = md_gpu_texture_storage(dst->color, 0);
    if (dst->surface) a.out_surface = md_gpu_texture_storage(dst->surface, 0);
    md_gpu_grid_t g = md_gpu_grid_for(v->r->k_export, v->width, v->height, 1);
    return md_gpu_launch(s, v->r->k_export, g, &a, sizeof(a));
}

void md_render_encode_id(uint32_t batch_index, uint64_t linear, uint32_t out_id[2]) {
    out_id[0] = (uint32_t)(linear & 0xFFFFFFFFu);
    out_id[1] = ((batch_index + 1u) << 20) | (uint32_t)((linear >> 32) & 0xFFFFFu);
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
    uint64_t linear = ((uint64_t)(w[1] & 0xFFFFFu) << 32) | w[0];
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

/* For md_render_bench: the compiled kernels, by index, with what the driver
   reports of each (registers, ...). False past the last. Not in the header. */
bool md_render_debug_kernel_info(md_render_ctx_t r, uint32_t i, const char** name, md_gpu_kernel_info_t* info) {
    if (!r) return false;
    const struct { const char* n; md_gpu_kernel_t k; } ks[] = {
        { "refit", r->k_refit }, { "cull_clusters", r->k_cull_clusters }, { "cull", r->k_cull }, { "finalize", r->k_finalize },
        { "raster_small", r->k_raster_small }, { "hiz", r->k_hiz }, { "resolve", r->k_resolve },
        { "mask_visible", r->k_mask_visible }, { "export_targets", r->k_export }, { "pick", r->k_pick },
        { "refit_heavy", r->k_refit_heavy }, { "cull_heavy", r->k_cull_heavy }, { "resolve_heavy", r->k_resolve_heavy }, { "resolve_poly", r->k_resolve_poly }, { "spline_order", r->k_spline_order },
    };
    if (i >= sizeof(ks) / sizeof(ks[0])) return false;
    if (name) *name = ks[i].n;
    if (info && !ks[i].k) { memset(info, 0, sizeof(*info)); return true; }   /* not made on this device: no stats */
    return info ? md_gpu_kernel_info(ks[i].k, info) : true;
}

/* The same for the raster pipelines: per stage, what the driver says */
bool md_render_debug_pipeline_info(md_render_ctx_t r, uint32_t i, const char** name, md_gpu_pipeline_info_t* info) {
    if (!r) return false;
    const struct { const char* n; md_gpu_pipeline_t p; } ps[] = {
        { "sphere", r->sphere_vis }, { "capsule", r->capsule_vis }, { "spline", r->spline_vis }, { "polygon", r->polygon_vis },
        { "sphere mask", r->sphere_mask }, { "capsule mask", r->capsule_mask }, { "spline mask", r->spline_mask }, { "polygon mask", r->polygon_mask },
    };
    if (i >= sizeof(ps) / sizeof(ps[0])) return false;
    if (name) *name = ps[i].n;
    if (info && !ps[i].p) { memset(info, 0, sizeof(*info)); return true; }
    return info ? md_gpu_pipeline_info(ps[i].p, info) : true;
}

