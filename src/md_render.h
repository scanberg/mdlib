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

        SPHERES     one point each
        CAPSULES    two points: a cylinder between them, round caps
                    (dashed, displaced); or flat ends, tapered (cones),
                    a radius in pixels (lines)
        POLYGONS    3..8 points around: a rounded slab
        SPLINES     chains of control points, swept stadium
        EXTERNAL    no primitives: a reserved id range for
                    something the caller renders itself

    Every primitive is raycast exactly; nothing is tessellated. Where a
    sphere and a capsule are the same surface (licorice joints), the sphere
    wins, whatever the order of the batches. Not yet: triangle meshes,
    capsules across a periodic boundary, transparency (below).

  * Points carry positions and, optionally, a VELOCITY. Through an optional
    index they also carry ATTRIBUTES (usually per atom): a colour or a palette
    index, two visibility bitsets and MARKS, a byte each of two state bits (an
    application's selection and highlight) and a tag (a highlight's group). A
    primitive is visible when all its points are, and has a state bit when all
    its points have it, with the largest of their tags.

  * A FRAME is the list of batches, a camera, a time stamp and three content
    counters. It is copied at the call. There is no scene graph and no object
    per molecule or per representation.

  * A VIEW owns the targets of one camera. Rasterisation writes depth and a
    64-bit id per pixel. The view then RESOLVES the surfaces (base colour,
    normal, motion) by intersecting again the one primitive that won each
    pixel, and makes MASKS of the two state bits (visible, and present at all
    including behind other surfaces), with their tags.

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
#include <core/md_bitset.h>

struct md_allocator_i;

#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
   Context
   ========================================================================= */

/* Pipelines and kernels shared by every view of one device. Not a global:
   tests render headless on a device of their own. */
typedef struct md_render* md_render_ctx_t;

typedef struct md_render_desc_t {
    struct md_allocator_i* alloc;   /* host allocations; NULL: the default heap */
    const char*            label;
} md_render_desc_t;

/* Synchronous and possibly slow (pipeline compilation): once, at start-up.
   NULL on failure (logged), including a device without graphics support. */
md_render_ctx_t md_render_create(md_gpu_device_t device, const md_render_desc_t* desc);

/* Every view must already be destroyed. Deferred like md_gpu's pipeline and
   kernel destruction. */
void md_render_destroy(md_render_ctx_t render);

/* =========================================================================
   Points
   =========================================================================

   The vertices of every kind of primitive. Positions (and velocities) are per
   point. Attributes are per attribute index: the point itself unless
   `attr_index` maps points onto another array (usually the atoms of a
   system). Every address but `position` is optional, and 0 has a defined
   meaning: no radius is 0, no colour is white, no palette index is entry 0,
   no visibility bitset is all visible, no marks are none set.

   The marks are the caller's memory, written by it: on unified memory straight
   from the CPU (md_render_marks_write fills a byte array from sets, for that or
   for an upload). md_render reads them when it makes the masks, so bump
   MD_RENDER_CONTENT_STATE when they change, and keep what a frame in flight
   reads unchanged until it is done (md_gpu_temp memory, or an upload). */
#define MD_RENDER_MARK_STATE0     0x01u
#define MD_RENDER_MARK_STATE1     0x02u
#define MD_RENDER_MARK_TAG_SHIFT  2
#define MD_RENDER_MARK_TAG_COUNT  8          /* tags 0 .. 7                         */
#define MD_RENDER_MARK_TAG(tag)   ((uint8_t)(((tag) & 7u) << MD_RENDER_MARK_TAG_SHIFT))

typedef struct md_render_points_t {
    uint32_t      count;

    md_gpu_addr_t position;        /* float3, `stride` bytes apart                    */
    uint32_t      stride;          /* multiple of 4; 0: 12                            */
    uint32_t      velocity_stride; /* multiple of 4; 0: 8                             */
    md_gpu_addr_t velocity;        /* packed half3 (see conventions); 0: at rest      */
    md_gpu_addr_t radius;          /* float per point (SPHERES), `radius_stride` bytes
                                      apart; 0: radius 0. With radius_table: uint32
                                      per point instead, an index into the table (past
                                      its count: radius 0). Read only by spheres with
                                      a radius_scale other than 0                     */
    uint32_t      radius_stride;   /* multiple of 4; 0: 4. Points stored as float4
                                      (x, y, z, radius or type): position with stride
                                      16 and radius = position + 12 with stride 16    */
    md_gpu_addr_t radius_table;    /* float per type: the radius of a point is
                                      radius_table[type]; 0: `radius` holds radii     */
    uint32_t      radius_table_count; /* entries                                    */

    uint32_t      attr_count;      /* entries of the attribute arrays; 0: count       */
    md_gpu_addr_t attr_index;      /* uint32 per point; 0: the point itself           */
    md_gpu_addr_t color;           /* uint32 RGBA8 per attribute, COLOR_POINTS        */
    md_gpu_addr_t palette_index;   /* uint16 per attribute, COLOR_PALETTE             */
    md_gpu_addr_t visible[2];      /* bitsets over the attributes, ANDed; 0: all set  */
    md_gpu_addr_t marks;           /* uint8 per attribute: bit 0 state 0, bit 1 state
                                      1, bits 2..4 the tag (MD_RENDER_MARK_*), 5..7
                                      zero; 0: none set                               */
} md_render_points_t;

/* The marks of `count` attributes into `dst`: the states from two sets (empty:
   none) and no tag. Then md_render_marks_tag for each set of a tag. */
void md_render_marks_write(uint8_t* dst, size_t count, md_bitset_t state0, md_bitset_t state1);
/* The tag of the attributes of `set` (the later call wins where sets overlap) */
void md_render_marks_tag(uint8_t* dst, size_t count, md_bitset_t set, uint32_t tag);

/* =========================================================================
   Batches
   ========================================================================= */

typedef enum md_render_kind_t {
    MD_RENDER_KIND_INVALID = 0,    /* a zero-initialised batch is an error */
    MD_RENDER_SPHERES,
    MD_RENDER_CAPSULES,
    MD_RENDER_POLYGONS,            /* rounded flat polygons: md_render_polygon_style_t */
    MD_RENDER_SPLINES,             /* tubes and ribbons along cubic splines: md_render_spline_style_t */
    MD_RENDER_EXTERNAL,            /* an id range for the caller's own renderer */
} md_render_kind_t;

typedef struct md_render_sphere_style_t {
    /* radius = point radius * radius_scale + radius_add; <= 0 is not drawn */
    float radius_scale;
    float radius_add;
} md_render_sphere_style_t;

typedef struct md_render_capsule_style_t {
    /* Radius of every capsule, times its own scale when `radius_scale` is
       given (a weak bond fading out; <= 0 hides it). Round caps: capsules
       join smoothly at a point. */
    float           radius;
    md_gpu_addr_t   radius_scale;  /* float per capsule, optional                      */

    /* Dashed (dash_length > 0): flat-ended cylinders of the radius between
       the spheres at its points (dash_ends: of radius point radius *
       radius_scale + radius_add, those a representation draws there),
       from where their rim meets each sphere, so that none is hidden or cut
       by one (md_gl's); zero dash_ends: from the points. Given by lengths,
       never by a count: as many dashes as come closest, at least two,
       stretched alike to fill the length. Nothing where not half a dash is
       left between the spheres. Only the hardware
       raster draws them (no compute raster for small ones). dash_length 0:
       solid. */
    float           dash_length;
    float           gap_length;
    md_render_sphere_style_t dash_ends;

    /* Displaced (offset != 0): each capsule moved off its points by
       `offset` times its radius, perpendicular to it, along its side: the
       line of `side` (float4 per capsule, xyz, in the points' space; 0: no
       side) weighed by its length against the line perpendicular to the
       capsule and the view, of weight 1/50, the two averaged as lines (by
       their doubled angles about the axis), the result towards `side`.
       So a side of length 1 holds, a short one turns smoothly towards
       facing the view, and none faces it. The thinner sticks of a multiple
       bond: batches over the same capsules and sides, the offset of either
       sign (md_gl's double and triple bonds).
       The ends of a displaced capsule pulled in: by its displacement times
       the miter of each end, side.w as two halves (packHalf2x16, the first
       end's low; 1 / tan of half the angle there, where a neighbour
       displaced alike meets it; 0: none), and by `trim` of its length;
       together at most four fifths of it. The dashed partners of
       delocalized bonds (md_gl's), an inner outline of their group. */
    float           offset;
    md_gpu_addr_t   side;
    float           trim;

    /* Shaped: cylinders and cones, and lines as wide on screen all along.
       `flat`: flat ends, not round (cylinders, truncated cones; the sharp
       edges of an arrow). `tapered`: the radius at the second point is
       end_radius (times the capsule's radius_scale, as the first's): a
       truncated cone, a cone where it is 0; with round ends the hull of
       the two balls. `radius_px` > 0: instead of the radii, at each end the
       world radius of that many pixels at its depth, so that a capsule is as
       wide on screen all along (a cell's edges; round ends close the joints
       of a polyline); a capsule that reaches past the near plane is cut
       short there. Not tapered. The radius_scale, the instance's scale and
       a fade's shrink multiply either. Not dashed or displaced. Only the
       hardware raster draws them (no compute raster for small ones). */
    bool            flat;
    bool            tapered;
    float           end_radius;
    float           radius_px;
} md_render_capsule_style_t;

/* SPLINES: a profile swept along a uniform cubic B-spline through control
   points, for cartoons, nucleic acid backbones, lipid tails and the like.
   The profile is a stadium (a 2D capsule): a segment of half length a
   along the profile's long axis, rounded by a radius r; a circle where a is
   0 (a tube). Formally, a part is the union of the balls of radius r around
   the swept segments.

   Control points: element e is control point e, at point index[e] & ~flags;
   the control points of a spline are consecutive elements, the first marked
   MD_RENDER_SPLINE_FIRST and the last MD_RENDER_SPLINE_LAST (both on one
   that is alone: not drawn). index is required.
   Parts: element e draws the part of its spline around its control point,
   from halfway to the one before to halfway to the one after (the first
   and the last only their inner half). A part takes the attribute of its
   control point: its colour, state bits, visibility and pick, so a hidden
   control point leaves a gap of its part. The spline passes near the
   control points (a B-spline smooths them) and through the first and the
   last (mirrored end points).
   Axes: the profile's long axis at each control point, any length (md_render
   normalises it and makes it orthogonal to the curve), interpolated along
   the spline as the positions are. Its sign is the caller's: the ribbon
   turns from one axis to the next the short way round only if they are
   within 90 degrees, so keep neighbours' signs consistent. Deciding that
   per frame from the axes alone flips pairs near 90 degrees apart back and
   forth with thermal noise; the decision needs memory across frames, which
   belongs above md_render (md_render_mol keeps it).
   Profile: a and r at the start and the end of each part, linear between:
   equal across a joint for a smooth spline (a coil tapering into a helix's
   ribbon over a part), unequal for a step (the base of an arrow). Ends are
   round: where a spline ends, next to a hidden part and at a step, a cap
   closes the part (the balls of its end's r around its end's segment).
   Through: a B-spline passes near its control points, smoothing them (off
   by a sixth of the second difference: 0.5 A on the C4' helix of DNA). With
   `through` each is moved first by one step of the inverse of that
   smoothing, from the two on either side (md_gl's tubes.through): the
   spline is off by the fourth difference over 36 instead, a sixteenth on
   DNA. For a spline something hangs off at its control points (the sugars
   of a nucleic acid off its C4').
   Fades: SHRINK scales a and r; THIN drops whole parts, which leaves their
   joints open: for splines, shrink.
   Not yet: flat caps. */
enum {
    MD_RENDER_SPLINE_FIRST = 1u << 31,
    MD_RENDER_SPLINE_LAST  = 1u << 30,
};

typedef struct md_render_spline_style_t {
    md_gpu_addr_t axis;            /* float4 per element (xyz), the long axis at
                                      its control point. 0: none (tubes)        */
    md_gpu_addr_t profile;         /* float4 per element: a where its part starts
                                      and where it ends, then r where it starts
                                      and where it ends. A part with r 0 at both
                                      ends is not drawn                         */
    float         radius;          /* without a profile: this r ...              */
    float         flat;            /* ... and this a (0: a tube)                */
    float         scale;           /* the profile's a and r times this; 0: 1     */
    bool          through;         /* through the control points, not near them */
} md_render_spline_style_t;

/* POLYGONS: a flat convex polygon made solid by a ball of `radius`: the
   Minkowski sum, a plate twice the radius thick with round edges (the rings
   of a cartoon, the sugars and bases of nucleic acids). With `plate` below 1
   the plate within the rims is thinner, its faces moved in to plate times
   the radius: a dish, the rims standing above it (md_gl's slabs.plate).
   Every polygon of a batch has n = `corners` corners (3 to 8): element e's
   are points n e .. n e + n - 1 in order around it (through index when it is
   given, as a capsule's two are). Rings of 5 and of 6 are two batches. Its plane is the corners' fit (Newell's normal, through
   their centroid); corners off it are projected onto it. The polygon is
   taken to be convex, as a ring is. A stick joining two polygons is a
   capsule of the same radius: the same surface as their edges.
   A polygon is visible, has a state bit, and survives a THIN fade when all
   its corners do; it takes its colour from its first corner, and a pick
   reports the corner nearest the hit as its attribute. SHRINK scales it
   about its centroid, radius and all.
   Raycast exactly: the plate's two faces, and a capsule along each edge.
   CIRCLES (inset > 0): in the polygon's stead the circle inside it, a
   torus of `radius` in its plane about its centroid (md_gl's circles of
   aromatic rings): its centre line `inset` inside the nearest edge (of
   radius the apothem less inset), none where that leaves less than twice
   `radius`. Coloured around it by the corners as a capsule is along it by
   its two (color.sharpness): between two corners by their angles about the
   centroid. `dashes` over each edge (0: solid), half of each edge's share
   of the circle, the gaps centred on the corners. `plate` does not apply.
   Raycast exactly (the quartic). */
typedef struct md_render_polygon_style_t {
    uint32_t corners;              /* 3 .. 8, the same for every polygon of the batch */
    float    radius;               /* half the thickness, the rounding of the edges;
                                      <= 0 not drawn                                 */
    float    plate;                /* the plate's half thickness within the rims, a
                                      fraction of the radius, up to 1; 0: 1, flat    */
    float    inset;                /* > 0: the circle inside it, so far from the edges */
    uint32_t dashes;               /* the circle's dashes over each edge, up to 16    */
} md_render_polygon_style_t;

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

/* Fading with distance. A band is a range of the projection's scale at an
   element (pixels per view-space unit at a sphere's centre or a capsule's
   midpoint), {start, end} with start > end > 0: where a transition begins,
   nearer, and where it ends, farther. Across `out` the batch fades out and
   beyond it is gone; across `in` it fades in and nearer it is absent (`in`
   nearer than `out`). {0, 0}: no such band. Through a band the transition is
   linear in view-space depth.
     SHRINK  the sizes (sphere and capsule radii) scale down to 0 across `out`
             and up from 0 across `in`: for continuous geometry, and a batch
             that fades to nothing.
     THIN    elements drop out one by one: an element is kept while a hash of
             its attribute (of both, for a capsule) is at least how far
             through `out` it is and below how far through `in`. A batch that
             fades in across the band another fades out across, with the same
             two numbers and over the same attributes, draws exactly the
             attributes the other leaves out: a handover without holes.
   Both may be set. Clusters wholly past a band are dropped by their bounds,
   unread. Which batches hand over to which is the caller's. The instance
   transform's scale does not enter; with an orthographic projection every
   element is equally far through a band. Elements dropped are counted in
   md_render_stats_t.faded.
   Keyed (SPHERES and CAPSULES): `key` gives each point another point of the
   batch's that stands for it (a group's centroid for its members), and that
   point decides its fade: its depth how far through the bands, its attribute
   what THIN hashes. A group fades as one, and a batch of the stand-ins
   themselves (unkeyed, over the same positions and attributes, with the same
   instances) fading the other way across the same band draws exactly the
   groups the keyed batch leaves out: a level of detail. A capsule is kept
   while both its ends are (each by its own key), and SHRINK scales it by the
   smaller of their factors. Clusters are dropped by the bounds of their
   visible elements: every key point must lie within them (a group's
   centroid does, when its whole group is in the cluster). */
enum {
    MD_RENDER_FADE_SHRINK = 1u << 0,
    MD_RENDER_FADE_THIN   = 1u << 1,
};

typedef struct md_render_fade_t {
    uint32_t      mode;            /* MD_RENDER_FADE_* or'ed; 0: no fading            */
    float         in[2];           /* {start, end}, start > end; {0, 0}: none          */
    float         out[2];
    md_gpu_addr_t key;             /* uint32 per point: the point that stands for it
                                      (an index into the same points); 0: itself  */
} md_render_fade_t;

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
           POLYGONS    n    0: points n i .. n i + n - 1 (n: style.polygons.corners)
           SPLINES     1    required: the control point and its flags
       A subset of the points is an index range or a visibility bitset. */
    md_gpu_addr_t      index;

    /* uint32 per value of `index` (or per index, without one): the point it
       stands for, so an index can hold something else the points are of
       (atoms, of points a regrouping keeps moving). 0: the index's values are
       points. Not SPLINES. */
    md_gpu_addr_t      index_map;

    /* uint32 per element, reported by picking instead of the element's index:
       the bond a capsule stands for. 0: the index. */
    md_gpu_addr_t      id;

    /* Bitsets (uint32 words, bit i of word i / 32) over the elements' ids
       (`id`, else their index), ANDed: an element whose bit is clear is not
       drawn, picked or marked, as one of an attribute that is not visible.
       Of what the elements stand for (the bonds of a frame, some of which
       broke). Part of the geometry (MD_RENDER_CONTENT_GEOMETRY). 0: all set.
       Not SPLINES. */
    md_gpu_addr_t      id_visible[2];

    /* Clusters: the first level of the cull. Cluster c is the elements
       cluster_offset[c] .. cluster_offset[c + 1] (uint32 [cluster_count + 1],
       in element order, as `index` orders them), less cluster_offset[0] (a
       batch over the tail of another's elements, its points moved to start
       there, can point into the middle of the other's offsets), best at most 64 and
       spatially compact: md_render_mol makes them from the topology for
       bonded structures and from a grid for solvent. 0: runs of 64
       consecutive elements, which is only as good as the element order.
       Their bounds are refitted when the geometry changes (not when only the
       camera moves), over the element centres, and padded by each batch's
       own style; batches with the same points, index, clusters, visibility
       and style share one refit. */
    md_gpu_addr_t      cluster_offset;
    uint32_t           cluster_count;

    union {
        md_render_sphere_style_t  spheres;
        md_render_capsule_style_t capsules;
        md_render_spline_style_t  splines;
        md_render_polygon_style_t polygons;
    } style;

    md_render_color_t  color;

    /* Instances: the batch drawn once per transform (periodic images,
       assemblies, systems side by side). Column-major float[16] each, copied
       at the call; rigid or uniformly scaled. NULL: once, as it is. */
    const float*       transform;
    uint32_t           transform_count;

    md_render_fade_t   fade;
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
    MD_RENDER_CONTENT_STATE,        /* the marks                                        */
    MD_RENDER_CONTENT_COUNT,
};

typedef enum md_render_jitter_t {
    MD_RENDER_JITTER_STATIC = 0,   /* jitter while nothing rasterised changes:
                                      progressive supersampling of a still image.
                                      Sample 0 is unjittered (pixel centres)       */
    MD_RENDER_JITTER_ALWAYS,       /* every draw rasterises a new sample: TAA       */
    MD_RENDER_JITTER_NONE,         /* one unjittered raster per change              */
    MD_RENDER_JITTER_GIVEN,        /* the caller's: frame.jitter_offset, a raster per
                                      change of it (an application's own TAA sequence,
                                      shared with what else it rasterises)         */
} md_render_jitter_t;

/* How a view rasterises. Every primitive is culled on the GPU each raster:
   by its visibility bits and radius, by the frustum, and by occlusion against
   a depth pyramid in two phases (the previous raster's pyramid first, then
   this raster's for what that found occluded, so the result is exact even
   when the camera or the geometry moved). Survivors go to the hardware
   raster, or, when they are a few pixels across, to a compute raster
   (64-bit atomics; the hardware raster when the device lacks them). */

typedef struct md_render_frame_t {
    uint32_t                 width, height;   /* of the targets, in pixels            */
    md_render_camera_t       camera;
    const md_render_batch_t* batches;         /* copied at the call; at most 4095     */
    uint32_t                 batch_count;
    md_render_jitter_t       jitter;
    uint32_t                 samples;         /* STATIC: rasters a still view makes
                                                 before it stops; 0 is taken as 1     */
    float                    jitter_offset[2];/* GIVEN: this raster's subpixel offset, in
                                                 pixels, x right and y down; motion is
                                                 without it, as without md_render's   */
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
md_render_view_t md_render_view_create(md_render_ctx_t render, const char* label);

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
                                  (also hidden), state 1 visible, state 1 anywhere;
                                  0 where not, else (8 + tag) / 15: the visible
                                  surface's tag, of anywhere the largest            */
} md_render_targets_t;

md_render_targets_t md_render_view_targets(md_render_view_t view);

/* ---- Export -----------------------------------------------------------------
   The targets into textures of the caller's, of formats every share takes
   (MD_GPU_TEX_SHARED, md_gpu_gl.h): for a compositor in another API to read
   the view, a fetch of each per pixel. Each optional (0: not written),
   STORAGE, of the targets' size:

       depth    R32_FLOAT     the depth target (reverse Z, 0: nothing there)
       color    RGBA8_UNORM   rgb: the colour target's; a: the mask target's
                              four channels as bits over 255 (1 state 0
                              visible, 2 state 0 anywhere, 4 state 1 visible,
                              8 state 1 anywhere), bits 4..6 the tag of
                              state 1 (visible, else anywhere), and 128 where
                              the visible surface's batch has marks
                              (points.marks): a surface of none has no state
                              to be without (a compositor dimming what is not
                              selected leaves it). Coverage is depth > 0
       surface  RGBA16_FLOAT  xy: the normal target's; zw: the motion target's

   Rows as the targets': the first is the top. As of the last draw; ordered in
   `stream` like any other work. */
typedef struct md_render_export_t {
    md_gpu_texture_t depth;
    md_gpu_texture_t color;
    md_gpu_texture_t surface;
} md_render_export_t;

bool md_render_view_export(md_gpu_stream_t stream, md_render_view_t view, const md_render_export_t* dst);

/* ---- Ids --------------------------------------------------------------------
       id.y  bits 31..20  1 + the batch's position in the frame (0: background)
             bits 19..0   linear index, high 20 bits
       id.x               linear index, low 32 bits
   linear = instance * count + element: 4095 batches per frame, 2^52 elements
   times instances per batch. An EXTERNAL batch's renderer writes ids made by
   md_render_encode_id with its own position in the frame. */

void md_render_encode_id(uint32_t batch_index, uint64_t linear, uint32_t out_id[2]);

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

/* ---- Statistics -------------------------------------------------------------
   What the cull did with the elements (times instances) of the last raster
   that has completed on the GPU. Asynchronous like picking: false until one
   has completed since the previous call. */

typedef struct md_render_stats_t {
    uint64_t elements;             /* submitted: count x instances over all batches */
    uint32_t hidden;               /* by visibility bits or a radius <= 0           */
    uint32_t outside;              /* the frustum                                   */
    uint32_t faded;                /* by their batch's fade                         */
    uint32_t occluded;             /* by the depth pyramid, after both phases       */
    uint32_t hardware[2];          /* rasterised in hardware: phase 1, phase 2      */
    uint32_t compute[2];           /* rasterised in compute:  phase 1, phase 2      */
    uint32_t occlusion_retested;   /* elements found occluded in phase 1, tested again */
    uint32_t clusters;             /* clusters x instances                          */
    uint32_t clusters_expanded[2]; /* whose elements were tested: phase 1, phase 2  */
    uint32_t clusters_occluded;    /* after both phases                             */
    bool     refitted;             /* the cluster bounds were refitted              */
    uint32_t raster;               /* which raster of the view these are from (1..) */
} md_render_stats_t;

bool md_render_view_stats(md_render_view_t view, md_render_stats_t* out_stats);

/* ---- Helpers ----------------------------------------------------------------- */

/* A velocity in the packed form md_render_points_t.velocity reads. */
void md_render_pack_velocity(const float v[3], uint32_t out[2]);

#ifdef __cplusplus
}
#endif

#endif /* MD_RENDER_H */
