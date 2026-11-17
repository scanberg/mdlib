/*
md_render_debug.h

Not part of the API: the switches md_render and md_render_mol keep for
measuring and testing (md_render_bench, the unit tests). None changes the
image but where it says so; most change only what it costs. Settings are
global unless they take a view.
*/

#ifndef MD_RENDER_DEBUG_H
#define MD_RENDER_DEBUG_H

#include <md_render.h>
#include <md_render_mol.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ---- md_render ---------------------------------------------------------------- */

/* A view's raster: without occlusion culling, without the compute raster,
   and the compute raster up to small_pixels across (0: 6). The next draw
   of the view rasterises. */
enum {
    MD_RENDER_DEBUG_NO_OCCLUSION     = 1u << 0,
    MD_RENDER_DEBUG_NO_COMPUTE_RASTER = 1u << 1,
};
void md_render_debug_view(md_render_view_t view, uint32_t options, float small_pixels);

/* Occlusion: phase 1's margin in radii, and the extent under which elements
   are not tested one by one (md_render.c says how they were chosen) */
void md_render_debug_occlusion(float margin_radii, float min_px);
/* Off: occlusion culling always runs, not by the view's measure */
void md_render_debug_occlusion_auto(bool on);
/* Pieces per box of a spline part in the hardware raster (above 16: the
   device's), and phase 1's spline boxes in this many buckets of depth,
   nearest first (0: the cull's order) */
void md_render_debug_spline_group(uint32_t pieces);
void md_render_debug_spline_order(uint32_t buckets);
/* 1: the spline and polygon fragments skip their hit and keep the proxy's
   depth (wrong pictures: what the vertex stage and the fill cost alone) */
void md_render_debug_raster_probe(uint32_t mode);
/* Called with each pass's name as it is recorded */
void md_render_debug_pass_hook(void (*hook)(void* user, const char* pass), void* user);
bool md_render_debug_kernel_info(md_render_ctx_t r, uint32_t i, const char** name, md_gpu_kernel_info_t* info);
bool md_render_debug_pipeline_info(md_render_ctx_t r, uint32_t i, const char** name, md_gpu_pipeline_info_t* info);

/* ---- md_render_mol ------------------------------------------------------------ */

void md_render_mol_debug_prefetch(uint32_t windows_ahead);
void md_render_mol_debug_write(uint32_t leave_out);
void md_render_mol_debug_frame_groups(uint32_t groups);
bool md_render_mol_debug_kernel_info(md_render_mol_ctx_t r, uint32_t i, const char** name, md_gpu_kernel_info_t* info);

#ifdef __cplusplus
}
#endif

#endif /* MD_RENDER_DEBUG_H */
