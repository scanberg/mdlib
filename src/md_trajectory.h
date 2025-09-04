#pragma once

#include <md_types.h>
#include <core/md_unit.h>

#include <stdint.h>
#include <stdbool.h>

enum {
	MD_TRAJECTORY_FLAG_NONE					= 0,
	MD_TRAJECTORY_FLAG_DISABLE_CACHE_WRITE	= 1,	
};

typedef uint32_t md_trajectory_flags_t;

struct md_trajectory_o;

typedef struct md_trajectory_header_t {
	size_t   num_frames;
	size_t   num_atoms;
	size_t   max_frame_data_size; // This represents the maximum size of any frame which is extracted using extract_frame_data.
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

typedef struct md_trajectory_i {
	struct md_trajectory_o* inst; // Opaque trajectory data

	bool (*get_header)(struct md_trajectory_o* inst, md_trajectory_header_t* header);

	// --- EASY MODE ---
    // Loads data for frame 'idx' into the supplied buffers. The buffers must be large enough to hold the data.
    // each field is optional and can be NULL if you don't need it.
	bool (*load_frame)(struct md_trajectory_o* inst, int64_t idx, md_trajectory_frame_header_t* header, float* x, float* y, float* z);

#if 0
	// --- ADVANCED MODE ---
	// Returns size in bytes of frame, data_ptr is optional and if supplied, the frame data will be written to it.
	size_t (*fetch_frame_data)(struct md_trajectory_o* inst, int64_t idx, void* data_ptr);

	// Decodes the raw frame data
	bool (*decode_frame_data)(struct md_trajectory_o* inst, const void* data_ptr, size_t data_size, md_trajectory_frame_header_t* header, float* x, float* y, float* z);
#endif
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

static inline size_t md_trajectory_max_frame_data_size(const md_trajectory_i* traj) {
	md_trajectory_header_t header;
	if (traj && traj->get_header && traj->get_header(traj->inst, &header)) {
		return header.max_frame_data_size;
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

static inline bool md_trajectory_load_frame(const md_trajectory_i* traj, int64_t idx, md_trajectory_frame_header_t* frame_header, float* x, float* y, float* z) {
	if (traj && traj->inst && traj->load_frame) {
		return traj->load_frame(traj->inst, idx, frame_header, x, y, z);
	}
	return false;
}

#if 0
static inline size_t md_trajectory_fetch_frame_data(const md_trajectory_i* traj, int64_t idx, void* data_ptr) {
	if (traj && traj->inst && traj->fetch_frame_data) {
		return traj->fetch_frame_data(traj->inst, idx, data_ptr);
	}
	return 0;
}

static inline bool md_trajectory_decode_frame_data(const md_trajectory_i* traj, const void* data_ptr, size_t data_size, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
	if (traj && traj->inst && traj->decode_frame_data) {
		return traj->decode_frame_data(traj->inst, data_ptr, data_size, header, x, y, z);
	}
	return false;
}
#endif
