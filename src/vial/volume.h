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
//   structure's own particles are not its target. Structures pair particle by particle: an array of sets alike (the
//   residues of one name), or set[k][w], k structures whose column j plays one part in each (smiles_atoms, mcs).
//   Kinds that are not alike, a tuple of such arrays, are each fitted onto a reference of their own, which is placed on
//   the first kind's by what the two have in common (md_mcs), or by poses=: the volume is the first kind's frame.
// - density(points, box=(lo, hi)): the points in a box of the run's frame, each taken in the image nearest the box's
//   centre. Without box=, over the unit cell, in fractional coordinates, so a cell that changes keeps its voxels; the
//   grid is that of frame 0's cell.

#include "interp.h"

#include <core/md_os.h>
#include <core/md_vec.h>

#define VIAL_VOLUME_MAX_VOXELS      ((uint64_t)512 * 512 * 512)
#define VIAL_VOLUME_SAMPLES         4.0     // per Å, by default
#define VIAL_VOLUME_ONE             16777216.0  // 2^24: a whole particle in the sum

// The extents of a grid of lengths len[3] at samples per Å: ceil(len * samples) each, rounded up to a multiple of 8.
// When that is more than VIAL_VOLUME_MAX_VOXELS, the samples per Å are lowered until it is not. Returns the samples per
// Å used; *coarser when they were lowered.
double vial_volume_grid(const double len[3], double samples, uint32_t dims[3], bool* coarser);

// Whether a built-in is one of the volumes, accumulated over the frames rather than computed at one
bool vial_volume_builtin(const char* name);

// One kind of structure of an sdf: count structures of n particles, column j playing one part in each, fitted by mass
// onto a reference of the kind, which R and t place on the first kind's
typedef struct vial_sdf_kind_t {
    uint32_t    count;
    uint32_t    n;
    uint32_t*   atoms;          // [count][n], in the order they pair
    uint32_t*   own;            // [count][n], each structure's ascending: not its own target
    md_float4_t*     ref_xyzw;       // [n]
    md_float3_t      ref_com;
    md_float3x3_t      R;              // the kind's reference about its centre to the first kind's about its centre
    md_float3_t      t;
} vial_sdf_kind_t;

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
    // An sdf: its kinds of structure, and the first kind's principal axes, the volume's frame
    struct vial_sdf_kind_t* kinds;
    uint32_t            num_kinds;
    uint32_t            num_structures; // over every kind, kind by kind
    uint32_t            max_n;
    md_float3x3_t              axes;
    const struct md_system_t* sys;
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
bool vial_volume_begin(vial_kctx_t* k, vial_volume_t* vol, const char* name, const vial_val_t* const* args, const md_str_t* names,
                       uint32_t n, const vial_frame_t* first, bool add, md_allocator_i* alloc);

// Adds the frame k->in->frame holds, from the op's arguments at that frame. Any thread, any order.
bool vial_volume_add(vial_kctx_t* k, vial_volume_t* vol, const vial_val_t* const* args, const md_str_t* names, uint32_t n);

// The volume over what was added: float[nx][ny][nz], smoothed by sigma, made in place of the sum, which out points at:
// it lives as long as vol. Once.
bool vial_volume_end(vial_kctx_t* k, vial_volume_t* vol, vial_val_t* out);

// A copy of vol set up alike, nothing added: for the same volume over other frames
bool vial_volume_like(vial_volume_t* out, const vial_volume_t* vol, md_allocator_i* alloc);

// Structure s of an sdf, over its kinds, kind by kind: its particles, ascending, n of them
const uint32_t* vial_volume_structure(const vial_volume_t* vol, uint32_t s, uint32_t* n);

// Structure s of an sdf at positions xyz in cell: the matrix from the run's frame to the volume's, column-major as
// md_float4x4_t is, the structure placed where it is drawn. False when it is not fitted at this frame, as when it spans too
// much of the cell to place.
bool vial_volume_pose(const vial_volume_t* vol, uint32_t s, const float* xyz, size_t num_atoms, const md_unitcell_t* cell, md_float4x4_t* out);

// The bulk density over the frames that had a cell, NaN when none had one
double vial_volume_bulk(const vial_volume_t* vol);

void vial_volume_free(vial_volume_t* vol);

// ### CORRESPONDENCE ###
// mcs(*structures): which particles of each structure play the parts of which of the first (md_mcs.h). Every
// structure is matched onto the first; the columns are the first's particles that every structure has a part for, in
// ascending order, and the result is set[k][w], a row per structure, one particle each. Structures alike (as many
// particles, of the same elements, or names, in order) share one mapping, found once. With xyz, the positions of
// every particle, of the mappings the bonds allow the one of least RMSD; without, the first. The check makes it
// without positions, evaluation with frame 0's, so the two differ at most by symmetric particles.
typedef struct vial_mcs_t {
    vial_sets_t*    rows;           // set[k][w], flat
    uint32_t        k, width;
    size_t          first_heavy;    // heavy atoms of the first structure
    uint32_t        least_at;       // the structure with the fewest pairs with the first (1-based; 0 when k is 1)
    size_t          least_common;   // its pairs, before the columns were cut to what all share
    size_t          least_heavy;    // the smaller heavy count of that pair
    bool            incomplete;
} vial_mcs_t;

bool vial_mcs(vial_mcs_t* out, const struct md_system_t* sys, const vial_sets_t* structures, bool by_name, const float* xyz, md_allocator_i* alloc);

// The built-in, from its arguments: every positional one a set or an array of sets, match= and least= by name.
// xyz as for vial_mcs, NULL for none.
bool vial_k_mcs(vial_kctx_t* k, const vial_val_t* const* args, const md_str_t* names, uint32_t n, const float* xyz, vial_val_t* out);
