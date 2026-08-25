#pragma once

#include <md_types.h>
#include <md_system.h>   // md_system_state_t is dereferenced below, a forward declaration is not enough
#include <core/md_unit.h>
#include <core/md_os.h>

#include <stdint.h>
#include <stdbool.h>

typedef enum {
	MD_TRAJECTORY_FLAG_NONE					= 0,
	MD_TRAJECTORY_FLAG_DISABLE_CACHE_WRITE	= 1,	
} md_trajectory_flags_t;

ENUM_FLAGS(md_trajectory_flags_t)

struct md_trajectory_o;
struct md_trajectory_reader_o;

// Coordinates, cell and frame ordinal are all handed over as a state, so that everything a frame
// yields travels together and cannot be paired with another frame's data.

typedef struct md_trajectory_header_t {
	size_t   num_frames;
	size_t   num_atoms;

	// Unit of frame_times. An EMPTY unit means the format could not tell us real time and the frame
	// times are ordinals standing in for it - monotonic and comparable, but not physical. Check it
	// before presenting a time to anyone.
	md_unit_t time_unit;

	// Physical time per frame, length num_frames. Never NULL.
	const double* frame_times;

	// Simulation step per frame, length num_frames, or NULL when the format does not record one.
	// This is the file's own notion of where a frame sits in the run - the LAMMPS TIMESTEP, the DCD
	// istart + i*nsavc - and it is NOT the frame ordinal. Indexed by ordinal, so
	// frame_steps[md_state_frame_nearest(state)] is the step a state was sampled at.
	const int64_t* frame_steps;
} md_trajectory_header_t;

// This is a common header we use for generated cache files
// Fileformats may have more fields, but this is the shared first part.
typedef struct md_trajectory_cache_header_t {
	uint64_t magic;
	uint64_t version;
	size_t num_bytes;
	size_t num_atoms;
	size_t num_frames;
	md_file_time_t last_modified; // modification timestamp of the source trajectory file
} md_trajectory_cache_header_t;

#ifdef __cplusplus
extern "C" {
#endif

typedef struct md_trajectory_reader_i {
	struct md_trajectory_reader_o* inst; // Opaque reader state
	void (*free)(struct md_trajectory_reader_i* self);

	// Everything a frame yields - coordinates, cell, atom count, ordinal - arrives on the state.
	// There is deliberately no second out param: a caller cannot pair one frame's metadata with
	// another frame's coordinates if there is only one thing to hold.
	bool (*load_frame)(struct md_trajectory_reader_o* inst, int64_t idx, md_system_state_t* state);
} md_trajectory_reader_i;

typedef struct md_trajectory_i {
	struct md_trajectory_o* inst; // Opaque trajectory data
	void (*free)(struct md_trajectory_i* self);

	// Retrieves the common trajectory header metadata.
	bool (*get_header)(struct md_trajectory_o* inst, md_trajectory_header_t* header);

	// Creates a reader with private I/O state suitable for reuse by a single worker/thread.
	// The trajectory object owns shared immutable metadata, while the reader owns transient state
	// such as file handles, scratch buffers and format specific decode state.
	bool (*init_reader)(md_trajectory_reader_i* reader, struct md_trajectory_o* inst);
} md_trajectory_i;

#ifdef __cplusplus
}
#endif

// Easy mode accessors
static inline size_t md_trajectory_num_frames(const md_trajectory_i* traj) {
	md_trajectory_header_t header;
	if (traj && traj->get_header && traj->get_header(traj->inst, &header)) {
		return header.num_frames;
	}
	return 0;
}

static inline size_t md_trajectory_num_atoms(const md_trajectory_i* traj) {
	md_trajectory_header_t header;
	if (traj && traj->get_header && traj->get_header(traj->inst, &header)) {
		return header.num_atoms;
	}
	return 0;
}

static inline md_unit_t md_trajectory_time_unit(const md_trajectory_i* traj) {
	md_trajectory_header_t header;
	if (traj && traj->get_header && traj->get_header(traj->inst, &header)) {
		return header.time_unit;
	}
	md_unit_t unit = {0};
	return unit;
}

static inline const double* md_trajectory_frame_times(const md_trajectory_i* traj) {
	md_trajectory_header_t header;
	if (traj && traj->get_header && traj->get_header(traj->inst, &header)) {
		return header.frame_times;
	}
	return NULL;
}

// NULL when the format does not record a simulation step per frame.
static inline const int64_t* md_trajectory_frame_steps(const md_trajectory_i* traj) {
	md_trajectory_header_t header;
	if (traj && traj->get_header && traj->get_header(traj->inst, &header)) {
		return header.frame_steps;
	}
	return NULL;
}

// Easy mode operations
static inline bool md_trajectory_get_header(const md_trajectory_i* traj, md_trajectory_header_t* header) {
	if (traj && traj->inst && traj->get_header) {
		return traj->get_header(traj->inst, header);
	}
	return false;
}

static inline void md_trajectory_free(md_trajectory_i* traj) {
	if (traj && traj->inst && traj->free) {
		traj->free(traj);
	}
}

// Reader mode operations
static inline bool md_trajectory_reader_init(md_trajectory_reader_i* reader, const md_trajectory_i* traj) {
	if (traj && traj->inst && traj->init_reader) {
		return traj->init_reader(reader, traj->inst);
	}
	return false;
}

static inline void md_trajectory_reader_free(md_trajectory_reader_i* reader) {
	if (reader && reader->free) {
		reader->free(reader);
	}
}

// Loads a frame into the supplied state: coordinates (when x/y/z are non NULL), the cell, the atom
// count, and the frame ordinal. To read only the metadata, hand over a state whose coordinate
// pointers are NULL - the readers skip the coordinates entirely.
//
// @NOTE: state->frame is stamped HERE rather than by the individual format readers. The ordinal is
// something the caller knows and the file does not, and the formats disagreed about it for years
// (some reported the simulation step, some the MODEL record number). Assigning it in one place
// makes that class of divergence unrepresentable rather than merely fixed. Do not push it back down
// into the readers.
static inline bool md_trajectory_reader_load_frame(md_trajectory_reader_i reader, int64_t idx, md_system_state_t* state) {
	if (reader.inst && reader.load_frame) {
		if (reader.load_frame(reader.inst, idx, state)) {
			if (state) {
				state->frame = (double)idx;
			}
			return true;
		}
	}
	return false;
}

static inline bool md_trajectory_load_frame(const md_trajectory_i* traj, int64_t idx, md_system_state_t* state) {
	bool result = false;
	if (traj && traj->inst && traj->init_reader) {
		md_trajectory_reader_i reader = {0};
		if (traj->init_reader(&reader, traj->inst)) {
			result = md_trajectory_reader_load_frame(reader, idx, state);
			md_trajectory_reader_free(&reader);
		}
	}
	return result;
}

// Physical time corresponding to a (possibly fractional) frame, in the trajectory's time_unit.
// This is the only place that decides what time means, including for formats which have none:
// md_trajectory_time_unit() reports an empty unit for those and the frame times fall back to
// ordinals, so the answer stays monotonic and comparable without pretending to be picoseconds.
//
// Returns -1.0 when the frame is absent (negative) or no frame times are available.
static inline double md_trajectory_time_at_frame(const md_trajectory_i* traj, double frame) {
	if (frame < 0.0) {
		return -1.0;
	}

	md_trajectory_header_t header;
	if (!traj || !traj->get_header || !traj->get_header(traj->inst, &header)) {
		return -1.0;
	}
	if (!header.frame_times || header.num_frames == 0) {
		return -1.0;
	}

	const double last = (double)(header.num_frames - 1);
	if (frame >= last) {
		return header.frame_times[header.num_frames - 1];
	}

	const size_t i = (size_t)frame;
	const double t = frame - (double)i;
	return header.frame_times[i] * (1.0 - t) + header.frame_times[i + 1] * t;
}

// The simulation step a (possibly fractional) frame sits at. A step is an integer count, so this
// does NOT interpolate - it reports the step of the nearest whole frame.
// Returns -1 when the frame is absent or the format records no steps.
static inline int64_t md_trajectory_step_at_frame(const md_trajectory_i* traj, double frame) {
	if (frame < 0.0) {
		return -1;
	}

	md_trajectory_header_t header;
	if (!traj || !traj->get_header || !traj->get_header(traj->inst, &header)) {
		return -1;
	}
	if (!header.frame_steps || header.num_frames == 0) {
		return -1;
	}

	size_t i = (size_t)(frame + 0.5);
	if (i >= header.num_frames) {
		i = header.num_frames - 1;
	}
	return header.frame_steps[i];
}
