#pragma once

#include <md_types.h>
#include <core/md_unit.h>

#include <stdint.h>
#include <stdbool.h>

typedef enum {
	MD_TRAJECTORY_FLAG_NONE					= 0,
	MD_TRAJECTORY_FLAG_DISABLE_CACHE_WRITE	= 1,	
} md_trajectory_flags_t;

ENUM_FLAGS(md_trajectory_flags_t)

struct md_trajectory_o;
struct md_trajectory_reader_o;

typedef struct md_trajectory_header_t {
	size_t   num_frames;
	size_t   num_atoms;
	md_unit_t time_unit;
    const double* frame_times;  // Array of length num_frames
} md_trajectory_header_t;

typedef struct md_trajectory_frame_header_t {
	size_t  num_atoms;
	int64_t index;
	double  timestamp;
	md_unitcell_t unitcell;
} md_trajectory_frame_header_t;

// This is a common header we use for generated cache files
// Fileformats may have more fields, but this is the shared first part.
typedef struct md_trajectory_cache_header_t {
	uint64_t magic;
	uint64_t version;
	size_t num_bytes;
	size_t num_atoms;
	size_t num_frames;
} md_trajectory_cache_header_t;

#ifdef __cplusplus
extern "C" {
#endif

typedef struct md_trajectory_reader_i {
	struct md_trajectory_reader_o* inst; // Opaque reader state
	void (*free)(struct md_trajectory_reader_i* self);

	bool (*load_frame)(struct md_trajectory_reader_o* inst, int64_t idx, md_trajectory_frame_header_t* header, float* x, float* y, float* z);
} md_trajectory_reader_i;

typedef struct md_trajectory_i {
	struct md_trajectory_o* inst; // Opaque trajectory data

	// Retrieves the common trajectory header metadata.
	bool (*get_header)(struct md_trajectory_o* inst, md_trajectory_header_t* header);

	// --- READER MODE ---
	// Creates a reader with private I/O state suitable for reuse by a single worker/thread.
	// The trajectory object owns shared immutable metadata, while the reader owns transient state
	// such as file handles, scratch buffers and format specific decode state.
	bool (*init_reader)(md_trajectory_reader_i* reader, struct md_trajectory_o* inst);

	// --- EASY MODE ---
    // Loads data for frame 'idx' into the supplied buffers. The buffers must be large enough to hold the data.
    // each field is optional and can be NULL if you don't need it.
	bool (*load_frame)(struct md_trajectory_o* inst, int64_t idx, md_trajectory_frame_header_t* header, float* x, float* y, float* z);
} md_trajectory_i;

typedef struct md_trajectory_loader_i {
	md_trajectory_i* (*create)(str_t filename, struct md_allocator_i* alloc, md_trajectory_flags_t flags);
	void (*destroy)(md_trajectory_i* traj);
} md_trajectory_loader_i;

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

// Easy mode operations
static inline bool md_trajectory_get_header(const md_trajectory_i* traj, md_trajectory_header_t* header) {
	if (traj && traj->inst && traj->get_header) {
		return traj->get_header(traj->inst, header);
	}
	return false;
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

static inline bool md_trajectory_reader_load_frame(md_trajectory_reader_i reader, int64_t idx, md_trajectory_frame_header_t* frame_header, float* x, float* y, float* z) {
	if (reader.inst && reader.load_frame) {
		return reader.load_frame(reader.inst, idx, frame_header, x, y, z);
	}
	return false;
}

static inline bool md_trajectory_load_frame(const md_trajectory_i* traj, int64_t idx, md_trajectory_frame_header_t* frame_header, float* x, float* y, float* z) {
	if (traj && traj->inst && traj->load_frame) {
		return traj->load_frame(traj->inst, idx, frame_header, x, y, z);
	}
	return false;
}
