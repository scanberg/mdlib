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

THE AMBIENT OCCLUSION is one of viamd's two scale-free screen-space
estimators, at half resolution, from a mip chain of the linear depth that
keeps one real pixel per texel (a rotated grid), with samples spread
log-uniformly between 2 pixels and a tenth of the view's height:
  QUALITY      horizon based with a visibility bitmask (ssao_bitmask.frag):
               2 slices of 8 steps a side; each sample occludes the arc from
               its front to its front pushed back by twice its distance
               beside the point; 32 sectors of equal shares of the
               cosine-weighted hemisphere. A thin occluder hides only the
               arc it covers. Visibility^(ao / 5)
  PERFORMANCE  obscurance (ssao.frag): 16 point samples, each judged at its
               own scale (its falloff radius 2.5 times its distance).
               1 - ao x the obscurance
Then a separable bilateral blur over two periods of their 4 x 4 pattern and
a joint bilateral upsample. The pattern is fixed: no history (a temporal
resolve after it accumulates what it would).

A transparent surface is a dielectric film of coverage a (its layer's
opacity), as viamd's isosurfaces were (isosurface.frag, isomesh.frag): of
what it is lit by it reflects a x diffuse + specular, so a clear surface
still shows its reflections, and lets through (1 - a)(1 - Fe), Fe the part
of the environment its specular lobe reflects (BASIC: a x its colour, and
1 - a). Front to back over what is behind them:
    C += T (a D + S);  T *= (1 - a)(1 - Fe)    (T: what is left, from 1)
and over the opaque surface of radiance R and coverage c (1, or 0 where
nothing was drawn) the result is (C + T c R, 1 - T (1 - c)).
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

/* ---- Volumes ------------------------------------------------------------------
   A scalar field on a 3D texture over a box, seen along each pixel's ray in
   front of the opaque surface, front to back with the layers: wherever a
   layer lies along the ray, what is in front of it is composited first. Each
   volume, cut to its clip box and to the kept side of its planes:
     EMISSION   what it emits and absorbs, through its transfer function
                (direct volume rendering), unshaded: the colour map's colour
     ISOSURFACES where the ray crosses one of its isovalues, a film shaded as
                a layer is (its colour's alpha the coverage), the optical
                density of the surfaces the ray is inside of absorbing the
                light in between. INSIDE isosurface i: the field on the far
                side of its value from 0 (f >= v for v >= 0, f <= v for v < 0)
   Marched at half a voxel of the finest volume the ray is in, all of them
   together (overlapping volumes are sampled at the same points); the
   emission's samples jittered by a pixel's hash (a temporal resolve after it
   averages them); a crossing found between two samples and refined on the
   trilinear field, so where it is does not depend on the sampling, but one
   the ray enters and leaves within a step is missed (as viamd's fast mode).
   An isosurface's colour is its iso_colors' times, when given, the colour
   volume's there, or the colour map's at the colour field's value there,
   with isolines of that field: lines where it takes given values, as wide in
   pixels however the surface is turned (|f - v| / |its change per pixel|, the
   change from the field's gradient and ray differentials), the persistent
   ones dark on a light rim, the transient one (hover) light on a dark rim,
   opaque paint on the surface. Its normal is the normals volume's direction
   there when given (and not too short), else the field's gradient.
   To come: the exact intersection (cell by cell). */

#define MD_COMPOSE_MAX_VOLUMES 4
#define MD_COMPOSE_MAX_ISO     8
#define MD_COMPOSE_MAX_PLANES  8
#define MD_COMPOSE_MAX_ISOLINES 8

typedef struct md_compose_volume_t {
    md_gpu_texture_t field;        /* 3D, R32_FLOAT or R16_FLOAT, SAMPLED: values at
                                      the voxels' centres, filtered trilinearly      */
    float    model[16];            /* its texture space [0, 1]^3 to the world
                                      (column-major)                                 */
    float    clip_min[3];          /* the box drawn, in texture space: 0 and 1 all */
    float    clip_max[3];
    uint32_t plane_count;
    float    planes[MD_COMPOSE_MAX_PLANES][4]; /* texture space, n . x + d >= 0 kept */

    /* Emission (none without a transfer function) */
    md_gpu_texture_t transfer;     /* 2D, N x 1, RGBA8_UNORM, SAMPLED: sRGB colour,
                                      the opacity of a slab 1/150 of the box thick  */
    float    transfer_range[2];    /* the field's values at its first and last texel */

    /* Isosurfaces */
    uint32_t iso_count;
    float    iso_values[MD_COMPOSE_MAX_ISO];
    float    iso_colors[MD_COMPOSE_MAX_ISO][4];   /* a: coverage                    */
    float    iso_densities[MD_COMPOSE_MAX_ISO];   /* optical density inside it:
                                                     extinction per world unit      */
    /* Optional, over the same box (each at its own resolution) */
    md_gpu_texture_t color_volume; /* 3D RGBA8_UNORM: times the isosurfaces' colours */
    md_gpu_texture_t color_field;  /* 3D R32_FLOAT or R16_FLOAT: through color_map
                                      (2D, N x 1, RGBA8_UNORM) over color_range; not
                                      with a colour volume                          */
    md_gpu_texture_t color_map;
    float    color_range[2];
    uint32_t isoline_count;        /* of the colour field                           */
    float    isoline_values[MD_COMPOSE_MAX_ISOLINES];
    bool     isoline_hover;        /* the transient line, at isoline_hover_value    */
    float    isoline_hover_value;
    float    isoline_width;        /* pixels                                        */
    md_gpu_texture_t normals;      /* 3D RGB10A2_UNORM or RGBA8_UNORM: the field's
                                      gradient's direction in world axes, (n + 1) / 2
                                      (md_gto_gpu_gradient_pack)                     */
} md_compose_volume_t;

/* ---- Slices --------------------------------------------------------------------
   A field shown on a rectangle in the world: where a pixel's ray meets it in
   front of the opaque surface, the transfer function's colour of the field's
   value there with the slice's opacity, composited in order along the ray
   with the layers and the volumes (viamd's planes display, dvr.frag's slice
   mode). Its edge: a line along the rectangle's border, edge_px pixels each
   way, on a dark halo a pixel wider, as wide in pixels wherever it is (the
   rectangle's size in pixels where the ray meets it, from the rays of the
   next pixels), in front of the slice. The field's value from a 3D texture
   through field_from_world, past the texture clamped to its edges: a
   volume's field (its texture space), or the slice's own image evaluated on
   the rectangle (W x H x 1, the rectangle's unit square as its x and y). */
#define MD_COMPOSE_MAX_SLICES 8

typedef struct md_compose_slice_t {
    float            rect[16];         /* the unit square x, y in [0, 1] (z 0) to the
                                          world (column-major)                       */
    md_gpu_texture_t field;            /* 3D, R32_FLOAT or R16_FLOAT, SAMPLED (null
                                          with opacity 0, as the transfer function) */
    float            field_from_world[16];
    md_gpu_texture_t transfer;         /* 2D, N x 1, RGBA8_UNORM: sRGB colour (its
                                          alpha unused)                              */
    float            transfer_range[2];
    float            opacity;          /* 0: its edge alone                          */
    float            edge[4];          /* sRGB colour, alpha: how much (0: none)     */
    float            edge_px;          /* half the line's width, pixels              */
} md_compose_slice_t;

/* ---- Shadows -------------------------------------------------------------------
   The key light's, on the opaque surface: from a shadow map, the depth of a view
   of the scene from the light (a depth-only md_render view, md_render_frame_t's
   depth_only, with an orthographic camera, reverse Z, whose +z is toward the
   light: the shading's light_dir, in this frame's view space), softened by the
   light's size (percentage-closer soft shadows: the blockers' mean distance
   from a search around the point, and the penumbra the light's angular radius
   makes over it, filtered with a disc of samples turned by each pixel and
   frame, which a temporal resolve averages), and refined in the frame's own
   depth by a short ray toward the light, `contact` of the map's texels long:
   what the map is too coarse for, where surfaces touch (within its biases,
   about two texels). The layers and the volumes are not shadowed. */
typedef struct md_compose_shadow_t {
    md_gpu_texture_t   depth;          /* that view's depth (md_render_targets_t)   */
    md_render_camera_t light;          /* its camera: orthographic                  */
    float              angular_radius; /* of the light, radians; 0: hard            */
    float              contact;        /* the ray's reach, texels; 0: none          */
} md_compose_shadow_t;

typedef enum md_compose_ao_mode_t {
    MD_COMPOSE_AO_QUALITY = 0,
    MD_COMPOSE_AO_PERFORMANCE,
} md_compose_ao_mode_t;

typedef struct md_compose_desc_t {
    md_render_camera_t   camera;   /* as the view was drawn (its projection: the
                                      direction each pixel looks along, and the
                                      depths' distances)                            */
    md_compose_shading_t shading;
    float                ao;       /* the occlusion's intensity (viamd's: 5); 0:
                                      none, not estimated                           */
    md_compose_ao_mode_t ao_mode;
    const md_compose_volume_t* volumes;    /* read at the call                      */
    uint32_t             volume_count;     /* at most MD_COMPOSE_MAX_VOLUMES        */
    const md_compose_slice_t*  slices;     /* read at the call                      */
    uint32_t             slice_count;      /* at most MD_COMPOSE_MAX_SLICES         */
    uint32_t             frame;            /* seeds the emission's jitter, and the
                                              shadows' samples                      */
    const md_compose_shadow_t* shadow;     /* null: none (read at the call)         */
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
