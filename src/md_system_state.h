#pragma once

#include <stdint.h>
#include <stdbool.h>
#include <md_types.h>
#include <md_unitcell.h>
#include <md_attributes.h>
#include <core/md_vec.h>

// The state of a system: the coordinates and the box of one frame (md_system_state_t).

// A snapshot of the geometric state of a system: where the atoms are and what box they are in.
//
// The FIELDS hold exactly what shares one interpolation contract - same type in and out, periodic
// boundary aware, handled as a unit by md_geom_interpolate_*. Nothing which interpolates differently
// becomes a field here, because a field in this struct is a promise that it does.
//
// 'attributes' is how a frame carries everything else. A TRR frame has velocities and forces beside
// its coordinates; they belong to that frame and to no other, and load_frame hands back one object
// precisely so a caller cannot pair one frame's velocities with another frame's positions. They are
// NOT part of the interpolation contract and nothing pretends otherwise: md_geom_interpolate_* takes
// raw coordinate arrays and never sees this table, so an interpolated state carries whatever its
// producer chose to put there and no quantity is silently blended.
//
// The table's allocator is the state's own, set by md_system_state_init, so a view state
// (alloc NULL) has no table - the same ownership rule the coordinates follow.
//
// What is computed from a frame's coordinates belongs there too: the backbone angles and the secondary structure of
// the frame (md_system_state_backbone_compute and the accessors beside it in md_util.h).
//
// Two presence bits, both self describing:
//   num_atoms == 0        -> no coordinates
//   unitcell.flags == 0   -> no cell
// num_atoms is used rather than testing xyz != NULL because md_array_ensure allocates capacity
// without setting size, so a non NULL xyz does not imply the coordinates were populated.
//
// xyz is packed, one md_float3_t of 12 bytes per atom - the layout every source and every consumer
// already has: the files, a run's atom/position attribute, the GPU buffers. Kernels that want x, y
// and z apart load four or eight atoms and de-interleave in registers (md_mm_load_xyz_packed_ps).
// The array is padded to a multiple of 16 atoms, zeroed past num_atoms, so such loads never need a
// scalar tail for reading.
//
// Ownership: alloc non-NULL means this state owns xyz and must be freed with
// md_system_state_free. alloc NULL means the state is a non owning view over coordinates somebody
// else owns - a scratch arena during script evaluation, a GPU mapped buffer during interpolation.
// Both forms are load bearing, and this is the only thing that distinguishes them.
//
// Set alloc before handing the state to a producer, exactly as md_system_t requires sys->alloc to
// be set before loading. The two then read symmetrically at the call site, and the state's lifetime
// is free to differ from the system's - a temp allocator for the state, a persistent one for the
// system, is a legitimate and useful combination.
// frame is the ordinal of the run frame the coordinates came from, as a continuous quantity: the
// integer part selects the frame, the fractional part is how far between that frame and the next the
// state has been interpolated. A NEGATIVE value means the state did not come from a run at all
// (a topology's own coordinates, a scratch buffer), which is why 0.0 cannot serve as that marker.
// Use md_state_has_frame / md_state_frame_floor / md_state_frame_nearest / md_state_frame_frac
// rather than reading the field: a raw cast of the absent value lands on frame 0, which is in range
// and plausible and therefore the kind of mistake that survives.
//
// CAVEAT, unlike the two presence bits above: -1 does not survive zero initialisation. A state
// built as {0} or with designated initialisers that omit frame reads as frame 0, not as absent.
// md_system_state_init stamps -1, so any state that went through it reads as absent until something
// writes a real frame; a hand rolled {0} does not. Only md_system_extract_frame and the
// interpolation which produces a state write a non negative value.
//
// @NOTE: physical time is deliberately NOT stored here. It is derivable from frame and the run's
// frame axis ("<run>/time"), and storing both would reintroduce the
// very thing this struct exists to prevent - two fields which must agree, with nothing enforcing it.
typedef struct md_system_state_t {
    size_t num_atoms;
    md_float3_t* xyz;
    md_unitcell_t unitcell;
    double frame;
    md_attributes_t attributes;   // per frame quantities beyond the interpolation contract above
    struct md_allocator_i* alloc;
} md_system_state_t;

#ifdef __cplusplus
extern "C" {
#endif

// State helpers

// A state which did not come from a trajectory carries a negative frame.
// @NOTE: NaN would pass this test, so nothing may ever write one into the field.
static inline bool md_state_has_frame(const md_system_state_t* state) {
    MD_ASSERT(state);
    return state->frame >= 0.0;
}

// The frame the state sits on or just after: the one to interpolate FORWARD from.
// @NOTE: a cast truncates toward zero, which is floor for the non negative values md_state_has_frame
// guarantees. Deliberately not floor() - keeping <math.h> out of this header matters, see below.
static inline int64_t md_state_frame_floor(const md_system_state_t* state) {
    MD_ASSERT(state);
    MD_ASSERT(md_state_has_frame(state));
    return (int64_t)state->frame;
}

// The frame the state most closely corresponds to. Use this where a single frame must be named -
// indexing a per frame array, for instance - not md_state_frame_floor, which would truncate 3.99 to 3.
static inline int64_t md_state_frame_nearest(const md_system_state_t* state) {
    MD_ASSERT(state);
    MD_ASSERT(md_state_has_frame(state));
    return (int64_t)(state->frame + 0.5);
}

// How far the state lies between md_state_frame_floor() and the frame after it, in [0,1).
static inline double md_state_frame_frac(const md_system_state_t* state) {
    MD_ASSERT(state);
    MD_ASSERT(md_state_has_frame(state));
    return state->frame - (double)(int64_t)state->frame;
}

static inline md_float3_t md_state_coord(const md_system_state_t* state, size_t atom_idx) {
    MD_ASSERT(state);
    if (atom_idx < state->num_atoms) {
        return state->xyz[atom_idx];
    }
    md_float3_t zero = {0};
    return zero;
}

// STATE
// A state owns its coordinate arrays and must be freed with the same allocator it was created with.

// True if the state carries coordinates. A state may legitimately carry a cell but no coordinates
// (a topology only format such as PSF), or coordinates but no cell (xyz without a cell).
static inline bool md_system_state_has_coords(const md_system_state_t* state) {
    return state && state->num_atoms > 0 && state->xyz;
}

static inline bool md_system_state_has_unitcell(const md_system_state_t* state) {
    return state && state->unitcell.flags != 0;
}

// Allocate coordinate storage for num_atoms using state->alloc, which must be set.
// Existing storage is freed first, so this doubles as the reset for a state being reused, and the
// padding up to the simd aligned capacity is zeroed. Returns false if the state has no allocator.
//
// CONTRACT for anything that produces a system and a state together (the format parsers):
// validate sys->alloc and state->alloc, then call md_system_reset(sys) and md_system_state_init(state, N)
// as a pair before touching either. Pass the exact atom count for N when it is known up front and
// write coordinates by index; pass 0 when atoms are filtered while parsing, then reserve with
// md_array_ensure(state->xyz, capacity, state->alloc) and push. Either way finish with
// state->num_atoms = sys->atom.count, and grow the coordinate array with state->alloc and never
// with sys->alloc - md_system_state_free releases them with state->alloc, and the two allocators
// are routinely different.
bool md_system_state_init(md_system_state_t* state, size_t num_atoms);

// Free the coordinate storage and zero the state, preserving the allocator so the state can be
// reused. A view (alloc == NULL) is simply zeroed. Takes no allocator argument by design: the state
// records the one it was allocated with, so it cannot be freed with the wrong one.
void md_system_state_free(md_system_state_t* state);

// Copy src into dst, reallocating dst as needed. dst->alloc must be set.
bool md_system_state_copy(md_system_state_t* dst, const md_system_state_t* src);

#ifdef __cplusplus
}
#endif
