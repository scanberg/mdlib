#pragma once

// A matrix of truths over the frames of a trajectory, rows × cols at each frame, kept as the frames each entry is true
// in, a bit per frame, for the entries true in some frame: a contact map over the trajectory, as the pairs ever in
// contact with when they are. Its size is the entries ever true times the frames, not rows × cols times the frames:
// 865 residues, 2,600 contacts a frame, perhaps three times as many pairs over a run, are 10 MB over 10,000 frames.
//
// Frames are added from any thread, each once, in any order, and read while they are: a reader sees the frames added
// so far, each whole. Over the frames an md_bitset_t gives (a run, a timeline period), an entry's fraction of them
// costs a popcount per 64 frames.

#include <core/md_bitset.h>

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

struct md_allocator_i;

#ifdef __cplusplus
extern "C" {
#endif

typedef struct md_pair_frames_t md_pair_frames_t;

// NULL when rows * cols is 0 or 2^32 or more, or num_frames is 0. Its memory from alloc, which must free.
md_pair_frames_t* md_pair_frames_create(uint32_t rows, uint32_t cols, uint32_t num_frames, struct md_allocator_i* alloc);
md_pair_frames_t* md_pair_frames_copy(const md_pair_frames_t* src, struct md_allocator_i* alloc);
void              md_pair_frames_free(md_pair_frames_t* pf);

uint32_t md_pair_frames_rows(const md_pair_frames_t* pf);
uint32_t md_pair_frames_cols(const md_pair_frames_t* pf);
uint32_t md_pair_frames_num_frames(const md_pair_frames_t* pf);

// Frame f: the entries true in it, i * cols + j each, in any order, duplicates allowed. False, and nothing added, when
// f is out of range or was added already. Entries out of range are left out.
bool md_pair_frames_add(md_pair_frames_t* pf, uint32_t frame, const uint32_t* entries, size_t count);

// Whether frame f was added; the frames added, with their words from alloc; the entries true in some frame added
bool        md_pair_frames_has(const md_pair_frames_t* pf, uint32_t frame);
md_bitset_t md_pair_frames_added(const md_pair_frames_t* pf, struct md_allocator_i* alloc);
size_t      md_pair_frames_count(const md_pair_frames_t* pf);

// The matrix at frame f into out[rows * cols], 1 where true and 0 elsewhere. False, and nothing written, when f was
// not added.
bool md_pair_frames_frame(uint8_t* out, const md_pair_frames_t* pf, uint32_t frame);

// Of the frames added among frames, the fraction each entry is true in, into out[rows * cols]. Returns how many frames
// that is; when none, out is NaN.
size_t md_pair_frames_occupancy(float* out, const md_pair_frames_t* pf, md_bitset_t frames);

// An entry true in some frame is kept in a slot, a number below md_pair_frames_count, for good. What is true in a part
// of the matrix, at a frame or over a period, is then a lookup in an index of the slots by entry, and a count of each:
// work as large as the part, which a contact map too large to go through whole at every frame shown can be drawn by,
// a tile at a time, on any thread.

// The entries true in some frame added, ascending, i * cols + j each, into entries, and the slot of each into slots:
// md_arrays from alloc, cleared first. Returns how many. It stays right while frames are added; what it lacks then
// are the entries first true in them, as many as md_pair_frames_count has grown by.
size_t md_pair_frames_index(uint32_t** entries, uint32_t** slots, const md_pair_frames_t* pf, struct md_allocator_i* alloc);

// Of the frames added among frames, how many each of slots[num] is true in, into counts[num]. Returns how many frames
// that is; a frame added meanwhile is in neither. Frames taken once from md_pair_frames_added make calls at different
// times count the same.
size_t md_pair_frames_slot_counts(uint32_t* counts, const md_pair_frames_t* pf, const uint32_t* slots, size_t num, md_bitset_t frames);

// The frames slot s is true in, with its words from alloc; empty when there is no such slot
md_bitset_t md_pair_frames_slot_frames(const md_pair_frames_t* pf, uint32_t slot, struct md_allocator_i* alloc);

// The frames entry (i, j) is true in, with its words from alloc; empty when it never is
md_bitset_t md_pair_frames_entry(const md_pair_frames_t* pf, uint32_t i, uint32_t j, struct md_allocator_i* alloc);

#ifdef __cplusplus
}
#endif
