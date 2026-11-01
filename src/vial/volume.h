#pragma once

// Volumes: sdf() and density(), a grid of a density accumulated over the frames (decided 1 November).
//
// A volume is never kept per frame. Each frame adds its particles to one running sum, spread over the eight voxels
// around each (cloud in cell), weighted by distance; the volume is the sum divided by what was added: frames, times
// structures for an sdf. The sum is of integers, each share rounded to 2^-24 of a particle: integers add up alike in
// any order, so a volume is the same bits whichever thread added which frame, as every value of the interpreter is. Its value is a density, per Å³ (or the weights' unit per Å³), absolute: published beside it
// is /bulk, the density of the same particles over the unit cell, so ρ/ρ_bulk is a division away.
//
// The grid follows from what it covers: samples per Å along each axis, each axis rounded up to a multiple of 8, and at
// most 512³ voxels, past which it is sampled coarser (W0407). Its voxels are in its own frame, written down as an
// origin, the corner of voxel [0][0][0], and three step vectors: the voxel [i][j][k] spans
// origin + (i..i+1) basis[0] + (j..j+1) basis[1] + (k..k+1) basis[2]. The value is float[nx][ny][nz], k fastest.
//
// - sdf(structures, target, extent): each structure is fitted, every frame, onto the reference (frame 0 of the first,
//   or ref=), with periodic boundaries, as v1 does; target's particles within the box of half extent `extent` about it
//   are added in the structure's frame. The volume's frame is the reference's centre of mass and principal axes. A
//   structure's own particles are not its target.
// - density(points, box=(lo, hi)): the points in a box of the run's frame, each taken in the image nearest the box's
//   centre. Without box=, over the unit cell, in fractional coordinates, so a cell that changes keeps its voxels; the
//   grid is that of frame 0's cell.

#include "interp.h"

#include <core/md_os.h>
#include <core/md_vec_math.h>

#define VIAL_VOLUME_MAX_VOXELS      ((uint64_t)512 * 512 * 512)
#define VIAL_VOLUME_SAMPLES         4.0     // per Å, by default
#define VIAL_VOLUME_ONE             16777216.0  // 2^24: a whole particle in the sum

// The extents of a grid of lengths len[3] at samples per Å: ceil(len * samples) each, rounded up to a multiple of 8.
// When that is more than VIAL_VOLUME_MAX_VOXELS, the samples per Å are lowered until it is not. Returns the samples per
// Å used; *coarser when they were lowered.
double vial_volume_grid(const double len[3], double samples, uint32_t dims[3], bool* coarser);

// Whether a built-in is one of the volumes, accumulated over the frames rather than computed at one
bool vial_volume_builtin(const char* name);

typedef enum vial_volume_kind_t {
    VIAL_VOLUME_SDF = 1,
    VIAL_VOLUME_BOX,
    VIAL_VOLUME_CELL,
} vial_volume_kind_t;

typedef struct vial_volume_t {
    uint8_t             kind;           // vial_volume_kind_t
    uint32_t            dims[3];
    double              origin[3];      // Å, in the volume's frame
    double              basis[3][3];    // Å: basis[a] is one step along index a
    double              samples;        // per Å, as used
    double              sigma;          // Å: a Gaussian the volume is smoothed with at the end; 0 for none
    md_unit_t           unit;           // of the value
    // An sdf: the structures, each of n particles, and the reference they are fitted onto, by mass
    const vial_sets_t*  structures;
    uint32_t            n;
    vec4_t*             ref_xyzw;
    vec3_t              ref_com;
    mat3_t              axes;           // the reference's frame to the volume's: its principal axes
    double              extent[3];      // the box's half extents
    // A box: its corner and lengths in the run's frame
    double              lo[3], len[3];
    // Over the cell: the cells added, summed, whose mean the published grid is laid out in
    double              cell_sum[3][3];
    // What is added, from any thread
    int64_t*            sum;            // [dims[0]][dims[1]][dims[2]], in VIAL_VOLUME_ONE per particle (or unit of weight)
    double              vref;           // Å³, a voxel's volume at frame 0: over a cell that changes, shares are scaled to it
    md_mutex_t          lock;
    bool                have_lock;
    uint64_t            added;          // frames, or frames times structures for an sdf
    uint64_t            frames;
    double              bulk;           // the target's density over the cell, summed over the frames with a cell
    uint64_t            bulk_frames;
    md_allocator_i*     alloc;
} vial_volume_t;

// Sets a volume up from its op's static arguments: its grid, and an sdf's reference, which is ref= or the first
// structure at frame `first` (the run's frame 0, which the caller reads; also the cell a density over the cell takes
// its grid from). With add, the sum frames are added to, zeroed, from the heap; without, only what places the volume, as
// drawing needs. False with an R diagnostic through k.
bool vial_volume_begin(vial_kctx_t* k, vial_volume_t* vol, const char* name, const vial_val_t* const* args, const str_t* names,
                       uint32_t n, const vial_frame_t* first, bool add, md_allocator_i* alloc);

// Adds the frame k->in->frame holds, from the op's arguments at that frame. Any thread, any order.
bool vial_volume_add(vial_kctx_t* k, vial_volume_t* vol, const vial_val_t* const* args, const str_t* names, uint32_t n);

// The volume over what was added: float[nx][ny][nz], smoothed by sigma, made in place of the sum, which out points at:
// it lives as long as vol. Once.
bool vial_volume_end(vial_kctx_t* k, vial_volume_t* vol, vial_val_t* out);

// A copy of vol set up alike, nothing added: for the same volume over other frames
bool vial_volume_like(vial_volume_t* out, const vial_volume_t* vol, md_allocator_i* alloc);

// Structure s of an sdf at positions xyz in cell: the matrix from the run's frame to the volume's, column-major as
// mat4_t is, the structure placed where it is drawn. False when it is not fitted at this frame, as when it spans too
// much of the cell to place.
bool vial_volume_pose(const vial_volume_t* vol, uint32_t s, const float* xyz, size_t num_atoms, const md_unitcell_t* cell, mat4_t* out);

// The bulk density over the frames that had a cell, NaN when none had one
double vial_volume_bulk(const vial_volume_t* vol);

void vial_volume_free(vial_volume_t* vol);
