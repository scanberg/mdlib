/*
md_render.h

Rendering of analytic primitives on md_gpu. Replaces md_gl (and retires
md_gfx). Molecular knowledge lives one level up, in md_render_mol.h; this
layer knows shapes, cameras and pixels.

    application        shading, AO, accumulation, tone mapping, UI
         |
    md_render_mol.h    molecules -> batches (bonds by kind, cartoons, ...)
         |
    md_render.h        batches of primitives -> view (visibility buffer,
         |             resolve, masks, picking, lazy rendering)
    md_gpu.h

The model, in full
------------------
  * A BATCH is a set of primitives of one KIND over a stream of POINTS that
    live in the caller's GPU memory, read in place:

        SPHERES     one point each                                [implemented]
        CAPSULES    two points: a cylinder between them           [implemented: round caps]
                    (cones, flat caps, dashes, periodic halves,
                    offsets: declared, not yet implemented)
        POLYGONS    2..7 points around: rounded slab or ring      [declared]
        SPLINES     chains of control points, swept stadium       [declared]
        TRIANGLES   meshes                                        [declared]
        EXTERNAL    no primitives: a reserved id range for
                    something the caller renders itself           [implemented]

    Every primitive but the triangle is raycast exactly; nothing is
    tessellated.

  * Points carry positions and, optionally, a VELOCITY. Through an optional
    index they also carry ATTRIBUTES (usually per atom): a colour or a palette
    index, two visibility bitsets and two state bitsets (an application's
    selection and highlight). A primitive is visible when all its points are,
    and has a state bit when all its points have it.

  * A FRAME is the list of batches, a camera, a time stamp and three content
    counters. It is copied at the call. There is no scene graph and no object
    per molecule or per representation.

  * A VIEW owns the targets of one camera. Rasterisation writes depth and a
    64-bit id per pixel. The view then RESOLVES the surfaces (base colour,
    normal, motion) by intersecting again the one primitive that won each
    pixel, and makes MASKS of the two state bits (visible, and present at all
    including behind other surfaces).

  * A view does the least work that keeps its targets correct:

        same batches, camera, size and counters      nothing (or the next
                                                     jittered sample, see below)
        only the STATE counter changed               masks
        only appearance changed                      resolve + masks, sample 0
        anything else                                raster + resolve + masks,
                                                     sample 0

    md_render cannot see writes into the caller's memory: the caller bumps the
    counter of the kind of content it changed. md_render_mol's helpers bump
    them for what they write.

    While nothing changes, successive draws rasterise along a fixed jitter
    sequence and report the sample index, so the caller can accumulate a
    supersampled image until it converges, after which nothing is drawn at
    all. Lazy rendering is how a view behaves, not a mode.

  * Motion is per view, from velocity and time. A view remembers the time,
    camera and instance transforms of its last raster; the displacement of a
    point is velocity * (frame.time - that time). So views drawn at different
    rates (insets, hidden views) are each correct without bookkeeping by the
    caller, and periodic images never streak: a velocity has no image.

  * md_render stops at surfaces. Shading, ambient occlusion, accumulation,
    tone mapping, volumes and overlays belong to the caller. Anything rendered
    by other means merges at the same targets, under an EXTERNAL batch's id.

  * Transparency (not implemented). The boundary is fixed so that the
    technique can change behind it: transparent batches resolve into N
    surface layers per pixel, front to back, which the caller shades and
    composites. Candidates: K-layer depth peeling into the visibility buffer
    (exact up to K, portable, no interlock); stochastic coverage under
    accumulation for heavy transparent geometry.

  * Ordering and memory follow md_gpu: every call that records GPU work takes
    the stream first, and nothing blocks. A view may resolve, mask or pick
    again from the arrays of the frame it last rasterised, without
    rasterising them again, so those arrays must stay valid and unchanged
    (other than as announced by the counters) until a later draw rasterises
    without them, or the view is destroyed. GPU temp scopes do not fit them.

Conventions
-----------
  World units are the caller's (mdlib: Ångström). Positions, radii and cells
  share them; velocities are world units per time unit of frame.time.

  View space looks down -Z (a right-handed look-at). Matrices are column-major
  float[16]. The projection is reverse Z: depth 1 at the near plane, 0 at the
  far plane, which may be infinite. Clip space follows md_gpu (+Y up, depth
  0..1). Perspective and orthographic projections both work.

  Colours are RGBA8 in a uint32_t, red in the low byte, colour channels
  sRGB-encoded. Bitsets are arrays of uint32_t: bit i is bit (i % 32) of
  word (i / 32).

  A velocity is three IEEE half floats packed in two uint32: x | y << 16,
  then z (high half unused). md_render_pack_velocity writes one.
*/

#ifndef MD_RENDER_H
#define MD_RENDER_H

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>

#include <core/md_gpu.h>

struct md_allocator_i;

#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
   Context
   ========================================================================= */

/* Pipelines and kernels shared by every view of one device. Not a global:
   tests render headless on a device of their own. */
typedef struct md_render* md_render_t;

typedef struct md_render_desc_t {
    struct md_allocator_i* alloc;   /* host allocations; NULL: the default heap */
    const char*            label;
} md_render_desc_t;

/* Synchronous and possibly slow (pipeline compilation): once, at start-up.
   NULL on failure (logged), including a device without graphics support. */
md_render_t md_render_create(md_gpu_device_t device, const md_render_desc_t* desc);

/* Every view must already be destroyed. Deferred like md_gpu's pipeline and
   kernel destruction. */
void md_render_destroy(md_render_t render);

/* =========================================================================
   Points
   =========================================================================

   The vertices of every kind of primitive. Positions (and velocities) are per
   point. Attributes are per attribute index: the point itself unless
   `attr_index` maps points onto another array (usually the atoms of a
   system). Every address but `position` is optional, and 0 has a defined
   meaning: no radius is 0, no colour is white, no palette index is entry 0,
   no visibility bitset is all visible, no state bitset is none set. */

typedef struct md_render_points_t {
    uint32_t      count;

    md_gpu_addr_t position;        /* float3, `stride` bytes apart                    */
    uint32_t      stride;          /* multiple of 4; 0: 12                            */
    uint32_t      velocity_stride; /* multiple of 4; 0: 8                             */
    md_gpu_addr_t velocity;        /* packed half3 (see conventions); 0: at rest      */
    md_gpu_addr_t radius;          /* float per point (SPHERES); 0: radius 0          */

    uint32_t      attr_count;      /* entries of the attribute arrays; 0: count       */
    uint32_t      _pad0;
    md_gpu_addr_t attr_index;      /* uint32 per point; 0: the point itself           */
    md_gpu_addr_t color;           /* uint32 RGBA8 per attribute, COLOR_POINTS        */
    md_gpu_addr_t palette_index;   /* uint16 per attribute, COLOR_PALETTE             */
    md_gpu_addr_t visible[2];      /* bitsets over the attributes, ANDed; 0: all set  */
    md_gpu_addr_t state[2];        /* bitsets over the attributes; 0: none set        */
} md_render_points_t;

/* =========================================================================
   Batches
   ========================================================================= */

typedef enum md_render_kind_t {
    MD_RENDER_KIND_INVALID = 0,    /* a zero-initialised batch is an error */
    MD_RENDER_SPHERES,
    MD_RENDER_CAPSULES,
    MD_RENDER_POLYGONS,            /* not yet implemented */
    MD_RENDER_SPLINES,             /* not yet implemented */
    MD_RENDER_TRIANGLES,           /* not yet implemented */
    MD_RENDER_EXTERNAL,            /* an id range for the caller's own renderer */
} md_render_kind_t;

typedef struct md_render_sphere_style_t {
    /* radius = point radius * radius_scale + radius_add; <= 0 is not drawn */
    float radius_scale;
    float radius_add;
} md_render_sphere_style_t;

typedef enum md_render_cap_t {
    MD_RENDER_CAP_ROUND = 0,       /* half sphere: capsules join smoothly at a point */
    MD_RENDER_CAP_FLAT,            /* not yet implemented                            */
} md_render_cap_t;

typedef enum md_render_offset_mode_t {
    MD_RENDER_OFFSET_WORLD = 0,    /* the offset vector as given                      */
    MD_RENDER_OFFSET_VIEW,         /* its length, perpendicular to both the capsule
                                      and the view direction (triple bonds)          */
} md_render_offset_mode_t;

typedef struct md_render_capsule_style_t {
    /* Radius of every capsule, times its own scale when `radius_scale` is
       given (a weak bond fading out; <= 0 hides it). */
    float           radius;
    md_render_cap_t cap[2];
    md_gpu_addr_t   radius_scale;  /* float per capsule, optional                      */

    /* Not yet implemented. Dashes are given by lengths, never by a count:
       as many flat-ended dashes as come closest (at least two), stretched
       alike. dash_length 0: solid. */
    float           dash_length;
    float           gap_length;

    /* Not yet implemented. A displacement of each capsule, in multiples of
       `radius` (float4 per capsule, xyz used): the thinner sticks of a double
       bond. With it the geometry no longer depends on the radius, so a
       generator runs once for all representations. */
    md_gpu_addr_t           offset;
    md_render_offset_mode_t offset_mode;

    /* Not yet implemented. Periodic capsules: each point joined to the
       nearest image of the other, drawn as two halves. All zero: none. */
    float           cell[9];
} md_render_capsule_style_t;

/* ---- Colour: only appearance; a change of it costs a resolve, not a raster */

typedef enum md_render_color_mode_t {
    MD_RENDER_COLOR_POINTS = 0,    /* points.color                                    */
    MD_RENDER_COLOR_UNIFORM,       /* `uniform`                                       */
    MD_RENDER_COLOR_PALETTE,       /* palette[points.palette_index]: one entry for a
                                      whole element or residue type                  */
} md_render_color_mode_t;

typedef struct md_render_color_t {
    md_render_color_mode_t mode;
    uint32_t      uniform;         /* RGBA8                                           */
    md_gpu_addr_t palette;         /* uint32 RGBA8 [palette_count]                    */
    uint32_t      palette_count;   /* indices past it take entry 0                    */
    float         sharpness;       /* CAPSULES: 0 blends the end colours smoothly, 1
                                      keeps each end's colour up to the middle       */
    float         tint[4];         /* rgb (linear), a: how much, 0 none               */
    float         desaturate;      /* 0 unchanged, 1 grey                             */
} md_render_color_t;

typedef struct md_render_batch_t {
    /* Identity across frames, chosen by the caller: non-zero, unique in a
       frame. Carries the previous instance transforms into motion, and comes
       back from picking. */
    uint64_t           key;
    md_render_kind_t   kind;

    /* Elements: spheres, capsules, ... 0 draws nothing. EXTERNAL: the number
       of ids the caller will write (informative). */
    uint32_t           count;
    md_render_points_t points;

    /* uint32 points per element:
           SPHERES     1    0: sphere i is point i
           CAPSULES    2    0: points 2i and 2i + 1
       A subset of the points is an index range or a visibility bitset. */
    md_gpu_addr_t      index;

    /* uint32 per element, reported by picking instead of the element's index:
       the bond a capsule stands for. 0: the index. */
    md_gpu_addr_t      id;

    union {
        md_render_sphere_style_t  spheres;
        md_render_capsule_style_t capsules;
    } style;

    md_render_color_t  color;

    /* Instances: the batch drawn once per transform (periodic images,
       assemblies, systems side by side). Column-major float[16] each, copied
       at the call; rigid or uniformly scaled. NULL: once, as it is. */
    const float*       transform;
    uint32_t           transform_count;
} md_render_batch_t;

/* =========================================================================
   Views
   ========================================================================= */

typedef struct md_render_view* md_render_view_t;

typedef struct md_render_camera_t {
    float view[16];                /* world to view                                   */
    float proj[16];                /* view to clip, reverse Z, without jitter         */
} md_render_camera_t;

/* Content counters: bump the counter of the kind of content changed in arrays
   a frame refers to, at unchanged addresses. Changes in the frame itself
   (addresses, counts, styles, colours, transforms, camera, size) are found by
   comparison. */
enum {
    MD_RENDER_CONTENT_GEOMETRY = 0, /* positions, velocities, radii, indices, visibility */
    MD_RENDER_CONTENT_APPEARANCE,   /* colours, palette indices, palettes               */
    MD_RENDER_CONTENT_STATE,        /* the state bitsets                                */
    MD_RENDER_CONTENT_COUNT,
};

typedef enum md_render_jitter_t {
    MD_RENDER_JITTER_STATIC = 0,   /* jitter while nothing rasterised changes:
                                      progressive supersampling of a still image.
                                      Sample 0 is unjittered (pixel centres)       */
    MD_RENDER_JITTER_ALWAYS,       /* every draw rasterises a new sample: TAA       */
    MD_RENDER_JITTER_NONE,         /* one unjittered raster per change              */
} md_render_jitter_t;

typedef struct md_render_frame_t {
    uint32_t                 width, height;   /* of the targets, in pixels            */
    md_render_camera_t       camera;
    const md_render_batch_t* batches;         /* copied at the call; at most 4095     */
    uint32_t                 batch_count;
    md_render_jitter_t       jitter;
    uint32_t                 samples;         /* STATIC: rasters a still view makes
                                                 before it stops; 0 is taken as 1     */
    double                   time;            /* of the positions, in the velocities'
                                                 time unit                           */
    uint64_t                 content[MD_RENDER_CONTENT_COUNT];
} md_render_frame_t;

typedef uint32_t md_render_work_t;
enum {
    MD_RENDER_WORK_RASTER  = 1u << 0,  /* depth, id                                  */
    MD_RENDER_WORK_RESOLVE = 1u << 1,  /* color, normal, motion                      */
    MD_RENDER_WORK_MASK    = 1u << 2,  /* mask                                       */
};

typedef struct md_render_status_t {
    md_render_work_t work;         /* 0: every target is as it was                    */
    uint32_t         sample;       /* index in the jitter sequence of what the targets
                                      hold; 0 after any change but of state          */
    float            jitter[2];    /* that raster's subpixel offset, in pixels         */
} md_render_status_t;

/* Records nothing: targets are made at the first draw, at its size. */
md_render_view_t md_render_view_create(md_render_t render, const char* label);

/* Deferred and non-blocking. */
void md_render_view_destroy(md_gpu_stream_t stream, md_render_view_t view);

/* Brings the targets up to date with the frame. On a GRAPHICS stream, outside
   any render pass. False on an invalid frame (logged): targets unchanged. */
bool md_render_view_draw(md_gpu_stream_t stream, md_render_view_t view,
                         const md_render_frame_t* frame, md_render_status_t* out_status);

/* ---- Targets ----------------------------------------------------------------
   Owned by the view, valid until a draw changes the size. All SAMPLED; id,
   color, normal and mask also STORAGE, so external renderers can merge. */

typedef struct md_render_targets_t {
    md_gpu_texture_t depth;    /* D32_FLOAT, reverse Z; 0 where nothing was drawn   */
    md_gpu_texture_t id;       /* RG32_UINT, see "Ids"                             */
    md_gpu_texture_t color;    /* RGBA8_UNORM, sRGB-encoded base colour; a: covered */
    md_gpu_texture_t normal;   /* RG16_FLOAT, view-space unit normal, octahedral    */
    md_gpu_texture_t motion;   /* RG16_FLOAT, pixels from where the surface was at
                                  the previous raster to where it is              */
    md_gpu_texture_t mask;     /* RGBA8_UNORM: state 0 visible, state 0 anywhere
                                  (also hidden), state 1 visible, state 1 anywhere */
} md_render_targets_t;

md_render_targets_t md_render_view_targets(md_render_view_t view);

/* ---- Ids --------------------------------------------------------------------
       id.y  bits 31..20  1 + the batch's position in the frame (0: background)
             bits 19..18  part (reserved: halves of periodic capsules, ...)
             bits 17..0   linear index, high 18 bits
       id.x               linear index, low 32 bits
   linear = instance * count + element: 4095 batches per frame, 2^50 elements
   times instances per batch. An EXTERNAL batch's renderer writes ids made by
   md_render_encode_id with its own position in the frame. */

void md_render_encode_id(uint32_t batch_index, uint32_t part, uint64_t linear, uint32_t out_id[2]);

/* ---- Picking ----------------------------------------------------------------
   Reads what the targets show after the last draw. Asynchronous: the result
   arrives once the GPU has got that far. */

typedef struct md_render_pick_t {
    uint64_t         key;          /* of the batch hit; 0: nothing there              */
    md_render_kind_t kind;
    uint32_t         batch_index;  /* in the frame last rasterised                    */
    uint32_t         instance;
    uint32_t         element;      /* in the batch                                    */
    uint32_t         id;           /* the batch's id of the element, else = element    */
    uint32_t         attr;         /* attribute index of the element's point nearest
                                      to the hit; UINT32_MAX for none (EXTERNAL)     */
    float            position[3];  /* world                                           */
    float            depth;
    uint32_t         x, y;         /* the pixel asked about                           */
} md_render_pick_t;

/* At most 4 requests in flight; a new one beyond that replaces the oldest. */
bool md_render_view_pick(md_gpu_stream_t stream, md_render_view_t view, uint32_t x, uint32_t y);

/* The most recent completed request not returned before. Does not block;
   completion is noticed through md_gpu_sync_is_complete. */
bool md_render_view_pick_result(md_render_view_t view, md_render_pick_t* out_pick);

/* ---- Helpers ----------------------------------------------------------------- */

/* A velocity in the packed form md_render_points_t.velocity reads. */
void md_render_pack_velocity(const float v[3], uint32_t out[2]);

#ifdef __cplusplus
}
#endif

#endif /* MD_RENDER_H */
