#pragma once

#include <stdint.h>
#include <stdbool.h>
#include <md_system.h>
#include <md_attributes.h>
#include <core/md_file.h>

// Runs: extracting a run's frames into a state, the run cache, and publishing runs and series.

#ifdef __cplusplus
extern "C" {
#endif

// EXTRACTION
//
// A state IS a snapshot of a run's temporal attributes at one frame, and these take it. What to take
// is said once, up front, and the context that comes back is then asked for frames:
//
//     const md_str_t paths[] = { MD_STR_LIT("atom/position"), MD_STR_LIT("unitcell") };
//     md_system_extract_t* ex = md_system_extract_begin(sys, run, paths, 2, alloc);
//     for (f = beg; f < end; ++f) md_system_extract_frame(ex, f, &state);
//     md_system_extract_end(ex);
//
// Saying it once is what lets the context be worth having. It resolves the paths a single time, and
// it keeps what a source is expensive to reopen - the files a trajectory streams from - open for as
// long as it lives. Opening a file is cheap on a local disk and not on the file servers of a
// cluster, where every open is a round trip to a metadata server that every rank shares.
//
// Paths are relative to the run ("<run>/<path>") and each lands in the state:
//
//   atom/position   into xyz, straight from the source: both are packed.
//   unitcell        into unitcell, a {F,3,3} box per frame with row i box vector i.
//   anything else   into out->attributes under the same relative path, as the value at that frame -
//                   a {F,N} c3 attribute arrives as {N} c3 - and without the temporal flag, since a
//                   snapshot has no frame axis. An attribute along an axis of its own (an energy file
//                   written more often than the coordinates, velocities written less often) is read
//                   at the row whose time matches the frame. A frame with no such row leaves it OUT
//                   of the state: its absence is the answer, and a state reused from frame to frame
//                   never shows another frame's value in its place.
//
// Every extract stamps out->frame. out->num_atoms must be 0 or the run's atom count, and a state
// asking for positions must have storage for that many; nothing is allocated for coordinates here.
//
// OWNERSHIP AND LIFETIME. The context belongs to the caller and to ONE thread at a time; a thread
// pool keeps one per thread. It holds attribute ids rather than pointers and looks each up again on
// every frame, so an attribute removed underneath it fails the extract instead of being read. It must
// still be ENDED before the run it reads is removed, since the files it holds are that run's.
typedef struct md_system_extract_t md_system_extract_t;

// NULL when the run, or any of the paths in it, does not exist. alloc is what the context lives in;
// the thread that begins it need not be the one that uses it, but only one may use it at a time.
md_system_extract_t* md_system_extract_begin(const md_system_t* sys, md_str_t run, const md_str_t paths[], size_t num_paths, struct md_allocator_i* alloc);
bool                 md_system_extract_frame(md_system_extract_t* ex, int64_t frame, md_system_state_t* out);
void                 md_system_extract_end(md_system_extract_t* ex);

// PUBLISHING A RUN
//
// What a trajectory format publishes is the same shape whatever the format (see RUNS):
//
//     <run>/time            F64 {F}        the frame axis; without a unit it holds frame ordinals
//     <run>/step            I64 {F}        the simulation step of each frame, when the file says
//     <run>/unitcell        F32 {F,3,3}    Angstrom, row i box vector i, zero where a frame has none
//     <run>/atom/position   F32 {F,N} c3   Angstrom, VIRTUAL: read from the file per frame
//     <run>/source/path     STR            the file the frames are read from
//     <run>/source/offset   I64 {F}        where each frame starts in it
//     <run>/source/size     I64 {F}        and how many bytes it has
//
// md_run_publish puts all of it in place from one description, or none of it: a publish that fails
// part way removes the run's prefix. What a format adds beyond this (velocities, how its frames are
// laid out) it publishes itself afterwards, under the same run, removing the prefix on failure too.
// Flags for a format's run publisher
typedef enum md_run_flag_t {
    MD_RUN_FLAG_NONE                = 0,
    MD_RUN_FLAG_DISABLE_CACHE_WRITE = 1,    // leave no '<file>.cache' index beside the file
} md_run_flag_t;

typedef uint32_t md_run_flags_t;

// INDEX CACHE. A format that has to scan its file to know where the frames are keeps what the scan
// learned beside it, as '<file>.cache', so the scan is paid once per file. The cache starts with
// this header and the format appends its own blocks after it.
//
// A cache is only as good as its match with the file: md_run_cache_open accepts it when the magic
// and version are the format's and the file's size and modification time, as the operating system
// reports them now, are the ones the cache was made from. Anything else - a run still being written,
// a file replaced by one of the same size - is a cache to make again.
typedef struct md_run_cache_header_t {
    uint64_t       magic;
    uint64_t       version;
    uint64_t       source_size;        // bytes
    md_file_time_t source_modified;    // the OS modification time, nanoseconds since the epoch
    uint64_t       num_atoms;
    uint64_t       num_frames;
} md_run_cache_header_t;

// Opens '<source_path>.cache' and reads its header. True, with the file left open just past the
// header for the format's own blocks, when the cache matches the source file as it is now.
bool md_run_cache_open(md_file_t* out_file, md_run_cache_header_t* out_header, md_str_t source_path, uint64_t magic, uint64_t version);

// Creates '<source_path>.cache' and writes the header, stamped with the source file's size and
// modification time as they were when it was scanned: take scanned with md_file_info_extract before
// the scan, so a file that grows during it leaves a cache that does not match it. True, with the
// file left open for the format's own blocks.
bool md_run_cache_create(md_file_t* out_file, md_str_t source_path, const md_file_info_t* scanned, uint64_t magic, uint64_t version, size_t num_atoms, size_t num_frames);

typedef struct md_run_desc_t {
    size_t          num_frames;
    size_t          num_atoms;

    const double*   time;               // num_frames values, required
    md_unit_t       time_unit;          // none when the file does not know time: time is then ordinals
    const int64_t*  step;               // num_frames values, or NULL when the file has none

    // The cell: RESIDENT from num_frames * 9 floats, or VIRTUAL from a provider (a format that
    // keeps it in the frames), or neither for a run without one.
    const float*                  unitcell;
    const md_attribute_virtual_t* unitcell_virt;

    md_str_t           source_path;        // copied
    const int64_t*  source_offset;      // num_frames values
    const int64_t*  source_size;        // num_frames values

    const md_attribute_virtual_t* position_virt;   // required
} md_run_desc_t;

bool md_run_publish(md_system_t* sys, md_str_t run, const md_run_desc_t* desc);

// "<run>/<leaf>" into buf. Empty when it does not fit.
md_str_t md_run_path(char* buf, size_t cap, md_str_t run, md_str_t leaf);

// For a provider of a run's attribute: the run it belongs to - its path minus "/<leaf>" - and the
// file and frame table the run reads from. False, and logged, when the attribute is not at <leaf>
// or the run has lost its source attributes.
typedef struct md_run_source_t {
    md_str_t          run;
    md_str_t          path;
    const int64_t* offset;
    const int64_t* size;
    size_t         num_frames;
} md_run_source_t;

bool md_run_source(md_run_source_t* out, const md_attributes_t* attributes, const md_attribute_t* attr, md_str_t leaf);

// SERIES ALONG A RUN
//
// Columns of numbers sampled along a run - an .xvg or .csv of per frame quantities - published as a
// group below it, "<run>/<group>", so that attr() and extraction read them like anything else in the
// run. Each column becomes "<run>/<group>/<name>", the name folded to lower case letters, digits and
// '_' (the original is the label); a name taken twice gets a suffix.
//
// With a time column the group has a frame axis of its own, "<run>/<group>/time", and every frame of
// the run must find its time in it, as an energy file's must; the columns are read at the matching
// row. Time without a unit is taken to be in the run's. Without a time column the rows ARE the run's
// frames, and there must be exactly as many.
//
// The group replaces one of the same name, all or nothing. source_path, when given, is published as
// "<run>/<group>/source" so a session can load it again.
typedef struct md_run_series_desc_t {
    md_str_t               group;          // below the run, e.g. "xvg/energy"
    size_t              num_rows;
    size_t              num_columns;
    const double*       time;           // num_rows values, or NULL: the rows are the run's frames
    md_unit_t           time_unit;
    const md_str_t*        names;          // num_columns
    const md_unit_t*    units;          // num_columns, or NULL for none
    const float* const* columns;        // num_columns arrays of num_rows values
    md_str_t               source_path;    // optional
} md_run_series_desc_t;

bool md_run_publish_series(md_system_t* sys, md_str_t run, const md_run_series_desc_t* desc);

#ifdef __cplusplus
}
#endif
