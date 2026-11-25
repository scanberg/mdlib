/*
md_post.h
What a frame's scene (md_compose.h) is made into to be shown, on md_gpu: the
display transform over the background, the temporal resolve (antialiasing,
and upscaling to the output's size) and the marks of the selection and the
highlight over the finished image, and the caller's overlay over that. Each
a call over textures of the caller's; nothing is kept between calls but what
the caller passes back (the history).

    md_compose.h        the scene: linear radiance, premultiplied, coverage
         |  md_post_display    exposure, curve, over the background, encoded
         |  md_post_resolve    TAA: this frame into the history, reprojected by
         |                     md_render's motion; at the output's size
         |  md_post_present    the marks over it, dithered: the image shown
         |  md_post_overlay    discs, segments, triangles over that
    the caller's present (an OpenGL blit, a swapchain)

Ported from viamd's post-processing (tonemap/display.frag, temporal.frag,
its marks pass), without the motion blur. The resolve runs on display
referred, encoded values (after the curve), as viamd's did: no HDR
fireflies in the history.
*/
#ifndef MD_POST_H
#define MD_POST_H

#include <stdint.h>
#include <stdbool.h>
#include <core/md_gpu.h>
#include <md_render.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct md_post* md_post_ctx_t;

md_post_ctx_t md_post_create(md_gpu_device_t device);
void          md_post_destroy(md_post_ctx_t post);

/* ---- Display transform ------------------------------------------------------
   x = exposure * L (1: reference white, the radiance of a white Lambertian
   surface under irradiance pi), shown as curve(x) / curve(W) with W three
   stops above reference white, sRGB-encoded and raised to 2.2 / gamma. The
   background lights nothing and is not curved: the scene is put over it after
   the curve by its coverage, so it shows as picked:
       display = a curve(x / a) / curve(W) + (1 - a) background             */
typedef enum md_post_curve_t {
    MD_POST_CURVE_NONE = 0,        /* clip to [0, 1] (basic shading, no tone mapping) */
    MD_POST_CURVE_NEUTRAL,         /* Khronos PBR Neutral                         */
    MD_POST_CURVE_ACES,            /* Stephen Hill's fit of the ACES RRT + ODT    */
    MD_POST_CURVE_HABLE,           /* John Hable's Uncharted 2 curve              */
} md_post_curve_t;

typedef struct md_post_display_t {
    md_gpu_texture_t scene;        /* md_compose's target                          */
    md_gpu_texture_t target;       /* RGBA16_FLOAT, STORAGE, of the scene's size:
                                      sRGB-encoded display values, a = 1           */
    md_post_curve_t  curve;
    float            exposure;     /* > 0                                          */
    float            gamma;        /* 2.2: the sRGB curve itself                   */
    float            background[3];/* sRGB, as picked                              */
} md_post_display_t;

bool md_post_display(md_gpu_stream_t stream, md_post_ctx_t post, const md_post_display_t* desc);

/* ---- Temporal resolve -------------------------------------------------------
   Playdead's TAA as viamd had it: the history reprojected by the motion of the
   nearest surface of the 3x3 (md_render's motion, which is the surface's that
   makes the most of a pixel), clipped toward the 3x3's colour box, blended by
   a feedback between feedback_min and feedback_max from the difference in
   luminance; where the 3x3 holds detail finer than a pixel the frame is taken
   from all nine samples (the prefilter). Upscaling (the history larger than the
   frame) blends each output pixel by how near a sample came to it, so a still
   image converges to the output's resolution. */
typedef struct md_post_resolve_t {
    md_render_view_t view;         /* the frame's depth and motion (its targets)   */
    md_render_camera_t camera;     /* as the view was drawn: the depths' distances */
    float            jitter[2];    /* this frame's subpixel offset, pixels, x right
                                      and y down (md_render_status_t.jitter)       */
    md_gpu_texture_t frame;        /* md_post_display's target                     */
    md_gpu_texture_t history_prev; /* the previous resolve's history, or 0: start
                                      over (the frame alone)                       */
    md_gpu_texture_t history;      /* RGBA16_FLOAT, STORAGE, SAMPLED: written, of
                                      the output's size (the frame's or larger)    */
    float            feedback_min; /* 0.80 */
    float            feedback_max; /* 0.95 */
} md_post_resolve_t;

bool md_post_resolve(md_gpu_stream_t stream, md_post_ctx_t post, const md_post_resolve_t* desc);

/* ---- Present ----------------------------------------------------------------
   The finished image (a history, or a display target) into the image shown,
   of the same size or larger (a frame upscaled without a resolve: bilinear),
   dithered to 8 bits, with the marks of the view over it: md_render's masks,
   filtered to the output's size. [0] the selection's (state 0), [1] the
   highlight's (state 1): a tint where in view, a fill where hidden behind
   something, an outline where present at all; each colour's alpha is how much
   of it. The highlight's tag t > 0 takes groups[t - 1] for all three.       */
typedef struct md_post_marks_t {
    float visible[2][4];
    float hidden[2][4];
    float outline[2][4];
    float groups[7][3];
    float highlight;               /* how much of the highlight shows: 0..1        */
    uint32_t outline_px;           /* the outline's width, output pixels (1..)     */
} md_post_marks_t;

typedef struct md_post_present_t {
    md_gpu_texture_t       image;  /* a history or a display target                */
    md_render_view_t       marks_view; /* its masks over the image; 0: none        */
    const md_post_marks_t* marks;
    md_gpu_texture_t       target; /* RGBA8_UNORM, STORAGE (and SHARED for a GL
                                      present): sRGB-encoded as the image is       */
    uint32_t               frame;  /* seeds the dither                             */
} md_post_present_t;

bool md_post_present(md_gpu_stream_t stream, md_post_ctx_t post, const md_post_present_t* desc);

/* ---- Overlay ----------------------------------------------------------------
   The caller's marks over the image shown, after md_post_present, into its
   target: annotations that the scene neither hides nor is hidden by, but where
   they are drawn (measurements, critical points, outlines). Drawn in order,
   each over what came before (the painter's), alpha-blended in the target's
   sRGB-encoded values, with the hardware raster at the target's size, of the
   camera unjittered.
     DISC      a flat disc at p0 facing the camera (in the view's x, y), of
               radius `size`, filled with `color` inside an outline `width`
               wide of `color2`. Rings (a highlight's) are larger discs before it.
     SEGMENT   from p0 to p1, `size` its half width, round ends, of `color`.
     TRIANGLE  p0, p1, p2, of `color`.
   Sizes in world units, or in pixels with MD_POST_OVERLAY_PIXELS (as wide
   wherever it is). Discs and segments are implicit: a quad each, whose
   fragments take their coverage from their distance to the disc's centre or
   the segment, in pixels, so their edges are smooth at any size; triangles
   are as rasterised. With MD_POST_OVERLAY_DEPTH an item is hidden where the
   scene is in front of it (`depth`: md_render's, of the same camera, at its
   own size; a disc at its centre's depth); otherwise it is over everything.
   Colours RGBA8, sRGB-encoded as md_render's, alpha how much. */
typedef enum md_post_overlay_kind_t {
    MD_POST_OVERLAY_DISC = 0,
    MD_POST_OVERLAY_SEGMENT,
    MD_POST_OVERLAY_TRIANGLE,
} md_post_overlay_kind_t;

enum {
    MD_POST_OVERLAY_PIXELS = 1u << 0,  /* size and width in pixels                   */
    MD_POST_OVERLAY_DEPTH  = 1u << 1,  /* hidden behind the scene                    */
};

typedef struct md_post_overlay_item_t {
    float    p0[3];
    uint32_t kind;                     /* md_post_overlay_kind_t                     */
    float    p1[3];
    uint32_t flags;
    float    p2[3];
    float    size;                     /* a disc's radius, a segment's half width   */
    uint32_t color;
    uint32_t color2;                   /* a disc's outline                           */
    float    width;                    /* the outline's, inside the radius; 0: none  */
    uint32_t _pad;
} md_post_overlay_item_t;

typedef struct md_post_overlay_t {
    md_gpu_texture_t   target;         /* RGBA8_UNORM, RENDER_TARGET: present's      */
    md_render_camera_t camera;         /* world to view and view to clip, unjittered */
    const md_post_overlay_item_t* items;   /* copied at the call                     */
    uint32_t           count;
    md_gpu_texture_t   depth;          /* the scene's (D32, reverse Z, SAMPLED) for
                                          MD_POST_OVERLAY_DEPTH; 0: none, drawn over */
} md_post_overlay_t;

bool md_post_overlay(md_gpu_stream_t stream, md_post_ctx_t post, const md_post_overlay_t* desc);

#ifdef __cplusplus
}
#endif

#endif /* MD_POST_H */
