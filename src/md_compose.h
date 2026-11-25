/*
md_compose.h
The transparent part of a frame on md_gpu: md_render's transparent surface
layers (md_render_targets_t.layers), shaded and composited front to back
into premultiplied linear radiance with coverage in alpha, for the caller to
put over its shaded opaque scene before tone mapping. Shading is md_render's
caller's; this is where it is done for what md_render leaves in layers.

    md_render.h       surfaces: the opaque targets, the transparent layers
         |
    md_compose.h      layers (and, to come, volumes between them) -> radiance

The shading models are those of viamd's deferred compose pass
(compose_deferred.frag), with the same parameters, so that transparent and
opaque surfaces of one colour under one light look alike:
  PBR    a simplified Cook-Torrance (GGX, Smith-Schlick, Schlick Fresnel) of a
         dielectric whose albedo is the colour; a key light, and an
         environment of two uniform hemispheres split at the view's horizon
         (sky above, ground below; view-space y): exact diffuse irradiance,
         and specular through Karis' analytic split-sum fit
  BASIC  the colours as given, darkened away from the light to `ambient` of
         themselves along a wrapped Lambert, (0.5 + 0.5 N.L)^2, with a
         Blinn-Phong highlight blended toward white
No ambient occlusion or shadows on transparent surfaces.

A layer of opacity a and shaded radiance L over what is behind it:
    C += T a L;  T *= 1 - a        (T: what is left of the view, from 1)
and the result is (C, 1 - T).
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

/* Kernels for one device */
md_compose_ctx_t md_compose_create(md_gpu_device_t device);
void             md_compose_destroy(md_compose_ctx_t compose);

typedef enum md_compose_shading_mode_t {
    MD_COMPOSE_SHADING_PBR = 0,
    MD_COMPOSE_SHADING_BASIC,
} md_compose_shading_mode_t;

/* View space throughout (the light's direction too) */
typedef struct md_compose_shading_t {
    md_compose_shading_mode_t mode;
    bool  linear_colors;           /* the layers' colours are linear reflectances,
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
                                      direction each pixel looks along)             */
    md_compose_shading_t shading;
    md_gpu_texture_t     target;   /* RGBA16_FLOAT, STORAGE, of the view's size:
                                      written, 0 where there is no layer            */
} md_compose_desc_t;

/* The layers of the frame `view` last drew, composited into desc->target. A
   view without layers (no transparent batch) writes 0 everywhere. False on a
   bad target (logged). */
bool md_compose(md_gpu_stream_t stream, md_compose_ctx_t compose, md_render_view_t view, const md_compose_desc_t* desc);

#ifdef __cplusplus
}
#endif

#endif /* MD_COMPOSE_H */
