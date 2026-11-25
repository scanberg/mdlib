/*
md_compose.h
The scene of a frame on md_gpu, from what md_render drew: the opaque surface
shaded under the ambient occlusion of its depth, and the transparent layers
(md_render_targets_t.layers) shaded and composited over it front to back,
into premultiplied linear radiance with the coverage in alpha (nothing
behind the scene: the background is the display's, md_post.h).

    md_render.h       surfaces: the opaque targets, the transparent layers
         |
    md_compose.h      surfaces (and, to come, volumes between them) -> radiance
         |
    md_post.h         radiance -> the image shown

The shading models are those of viamd's deferred compose pass
(compose_deferred.frag), with the same parameters, so that a transparent and
an opaque surface of one colour under one light look alike:
  PBR    a simplified Cook-Torrance (GGX, Smith-Schlick, Schlick Fresnel) of a
         dielectric whose albedo is the colour; a key light, and an
         environment of two uniform hemispheres split at the view's horizon
         (sky above, ground below; view-space y): exact diffuse irradiance,
         and specular through Karis' analytic split-sum fit. The occlusion
         takes the environment's diffuse light, and its specular light
         through Lagarde's specular occlusion
  BASIC  the colours as given, darkened away from the light to `ambient` of
         themselves along a wrapped Lambert, (0.5 + 0.5 N.L)^2, with a
         Blinn-Phong highlight blended toward white; the occlusion darkens
         all of it (the depth cue it is there)
No occlusion on transparent surfaces, and no shadows anywhere (yet).

THE AMBIENT OCCLUSION is viamd's scale-free screen-space obscurance (its
"performance" estimator): at half resolution, 16 samples spread
log-uniformly between 2 pixels and a tenth of the view's height, each judged
at its own scale (its falloff radius 2.5 times its distance), from a mip
chain of the linear depth that keeps one real pixel per texel (a rotated
grid); a separable bilateral blur over two periods of its 4 x 4 pattern; a
joint bilateral upsample. The pattern is fixed: no history (a temporal
resolve after it accumulates what it would).

A layer of opacity a and shaded radiance L over what is behind it:
    C += T a L;  T *= 1 - a        (T: what is left of the view, from 1)
and over the opaque surface of radiance S and coverage c (1, or 0 where
nothing was drawn) the result is (C + T c S, 1 - T (1 - c)).
*/
#ifndef MD_COMPOSE_H
#define MD_COMPOSE_H

#include <stdint.h>
#include <stdbool.h>
#include <core/md_gpu.h>
#include <md_render.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct md_compose* md_compose_ctx_t;

/* Kernels for one device. The occlusion's scratch is temp memory of the
   stream, for the call. */
md_compose_ctx_t md_compose_create(md_gpu_device_t device);
void             md_compose_destroy(md_compose_ctx_t compose);

typedef enum md_compose_shading_mode_t {
    MD_COMPOSE_SHADING_PBR = 0,
    MD_COMPOSE_SHADING_BASIC,
} md_compose_shading_mode_t;

/* View space throughout (the light's direction too) */
typedef struct md_compose_shading_t {
    md_compose_shading_mode_t mode;
    bool  linear_colors;           /* the surfaces' colours are linear reflectances,
                                      not sRGB                                       */
    float light_dir[3];            /* toward the key light (normalised here)        */
    float light_strength;          /* of the key light (both modes)                 */
    /* PBR */
    float dir_radiance[3];         /* the key light: irradiance at normal incidence */
    float sky_radiance[3];         /* the environment above the view's horizon      */
    float ground_radiance[3];      /* and below                                     */
    float roughness;               /* clamped to [0.04, 1]                          */
    float f0;                      /* 0.04: plastic                                 */
    /* BASIC */
    float ambient;
    float specular;
} md_compose_shading_t;

typedef struct md_compose_desc_t {
    md_render_camera_t   camera;   /* as the view was drawn (its projection: the
                                      direction each pixel looks along, and the
                                      depths' distances)                            */
    md_compose_shading_t shading;
    float                ao;       /* the occlusion's intensity: the visibility is
                                      1 - ao x the obscurance (viamd's: 5); 0: none,
                                      not estimated                                 */
    md_gpu_texture_t     target;   /* RGBA16_FLOAT, STORAGE, of the view's size:
                                      the scene's radiance premultiplied by its
                                      coverage (alpha), 0 where nothing was drawn   */
} md_compose_desc_t;

/* The scene of the frame `view` last drew, into desc->target. False on a bad
   target (logged). */
bool md_compose(md_gpu_stream_t stream, md_compose_ctx_t compose, md_render_view_t view, const md_compose_desc_t* desc);

#ifdef __cplusplus
}
#endif

#endif /* MD_COMPOSE_H */
