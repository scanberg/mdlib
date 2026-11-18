#pragma once

// PIECES: a long GPU evaluation as a sequence of bounded submissions.
//
// A submission that runs for seconds loses the device to a watchdog (TDR on Windows, 2 s; macOS ends command
// buffers that keep the display from the GPU), and before that keeps every other context's work - the
// rendering, in OpenGL - waiting behind it. No API gives a portable guarantee better than this: the driver
// arbitrates between contexts at submission boundaries, and at best at dispatch boundaries (mid-dispatch compute
// preemption exists on some GPUs and OSes, not on all of the ones mdlib runs on). A separate compute queue gives
// overlap, not isolation: a dispatch that fills the GPU with long-lived groups starves the display until they
// retire. So a long evaluation goes as PIECES: each a range of its work, each its own submission, and only a
// little of it issued ahead of what the GPU has finished.
//
// This is the part every such evaluation shares, the kernels' own code being the rest: how much the next piece
// takes, how much is issued ahead, the submission and its marker, how far the work is, and the piece's GPU time.
// A piece is sized to a target in milliseconds of GPU time (piece_ms), from the rate the pieces before it ran at
// as the device measured them (md_gpu_sync_gpu_seconds), so the same code makes pieces of a few ms on a GTX 1060,
// an RTX 3080, an Intel iGPU and an M4 without a cost model of its kernels. Before the first piece is measured,
// and on a queue that cannot time, it is the caller's estimate (first_units).
//
// The work is counted in UNITS, whatever the caller's are: points, groups, bricks - or better, a work estimate,
// which makes the measured rate steadier where the cost of a point varies (empty space against the nuclei).
// A piece is a contiguous range [issued, issued + n) of them.
//
// Embedded in a job; a poll once a frame:
//
//     md_gpu_pieces_poll(&p);                                   // what has passed, measured
//     for (uint64_t n; (n = md_gpu_pieces_next(&p)) > 0;) {     // room for another piece: n units
//         issue units [p.issued, p.issued + n) into p.stream;
//         md_gpu_pieces_end(&p, n);                             // submitted, marked
//     }
//     if (md_gpu_pieces_done(&p)) ...                           // every piece has passed
//
// A caller that must have it all at once (a blocking launch) sets ahead_ms < 0: every piece is issued, each its
// own submission still, and md_gpu_pieces_wait waits for them.
//
// One thread at a time: the stream's (md_gpu.h). Nothing here allocates.

#include <stdbool.h>
#include <stdint.h>

#include <core/md_gpu.h>

#ifdef __cplusplus
extern "C" {
#endif

// For measurements: bumped with every change to how pieces are sized or issued
#define MD_GPU_PIECES_REVISION "pieces r1 (4 ms a piece by the measured GPU rate, 24 ms ahead, doubled up to 96 ms when the GPU runs dry)"

#define MD_GPU_PIECES_PIECE_MS 4.0
#define MD_GPU_PIECES_AHEAD_MS 24.0
#define MD_GPU_PIECES_MARKS    64

typedef struct md_gpu_pieces_policy_t {
    // The GPU time a piece aims at. 0: MD_GPU_PIECES_PIECE_MS. The longest the display waits behind one piece:
    // a few ms costs a frame little, and a submission costs tens of microseconds.
    double piece_ms;
    // The GPU time issued ahead of what has passed. 0: MD_GPU_PIECES_AHEAD_MS. Enough that the GPU does not run
    // dry between two polls (a frame apart); each poll that finds it did doubles it, up to ahead_ms_max. Less
    // ahead: the display gets the GPU sooner, the evaluation takes longer. < 0: no bound, everything is issued
    // at once (in pieces still).
    double ahead_ms;
    double ahead_ms_max;    // 0: 4 x ahead_ms
} md_gpu_pieces_policy_t;

typedef struct md_gpu_pieces_mark_t {
    md_gpu_sync_t sync;
    uint64_t      end;      // the units issued once this piece is
    uint64_t      units;    // of this piece
} md_gpu_pieces_mark_t;

typedef struct md_gpu_pieces_t {
    md_gpu_stream_t stream;
    uint64_t total;         // units in all
    uint64_t issued;        // units issued: the next piece starts here
    uint64_t done;          // units whose pieces have passed

    // sizing
    uint64_t first_units;   // a piece before any is measured, or on a queue that cannot time
    uint64_t min_units;     // granularity: every piece but the last a multiple of it
    uint64_t max_units;     // the most a piece may take (the caller's scratch)
    uint64_t last_units;    // the last piece's
    double   rate;          // units per ms of GPU time, measured; 0 until then
    double   piece_ms, ahead_ms, ahead_ms_max;
    bool     untimed;       // a passed piece had no GPU time: the queue cannot time, sized by first_units

    // the pieces in flight, oldest at head
    md_gpu_pieces_mark_t mark[MD_GPU_PIECES_MARKS];
    uint32_t head, tail;

    // the pieces passed whose GPU time was not there yet (Metal's comes a moment later): asked again
    md_gpu_pieces_mark_t late[8];
    uint8_t  late_polls[8];
    uint32_t num_late;

    // what it took
    uint32_t num_pieces;
    uint32_t num_measured;
    double   gpu_ms;        // the measured pieces' GPU time
    uint64_t measured_units;
    double   max_piece_ms;  // the longest piece measured
    uint32_t num_dry;       // polls that found the GPU had run dry
} md_gpu_pieces_t;

// total: the units of the evaluation. first_units: the caller's estimate of a piece of piece_ms (0: total / 64).
// min_units, max_units: a piece's granularity and its largest (0: 1 and total). policy may be NULL.
void md_gpu_pieces_init(md_gpu_pieces_t* p, md_gpu_stream_t stream, uint64_t total, uint64_t first_units, uint64_t min_units, uint64_t max_units,
                        const md_gpu_pieces_policy_t* policy);

// The units the next piece may take now: 0 when all is issued, or enough is ahead of the GPU (or, bounded, the
// marker ring is full). Otherwise at least min_units (or what is left), at most max_units.
uint64_t md_gpu_pieces_next(md_gpu_pieces_t* p);

// The caller has issued the units [issued, issued + units) into the stream: submitted (md_gpu_stream_record),
// marked, counted. False when nothing could be submitted.
bool md_gpu_pieces_end(md_gpu_pieces_t* p, uint64_t units);

// Takes the pieces that have passed, and their GPU times for the rate. True when every piece is done. Never waits.
bool md_gpu_pieces_poll(md_gpu_pieces_t* p);

// Waits for the oldest piece in flight (none: returns at once)
void md_gpu_pieces_wait(md_gpu_pieces_t* p);

static inline bool md_gpu_pieces_all_issued(const md_gpu_pieces_t* p) { return p->issued >= p->total; }
static inline bool md_gpu_pieces_done(const md_gpu_pieces_t* p) { return p->done >= p->total; }

// How far the evaluation is, 0 to 1: the units of the pieces that have passed
static inline float md_gpu_pieces_progress(const md_gpu_pieces_t* p) {
    return p->total ? (float)((double)p->done / (double)p->total) : 1.0f;
}

// One line of what it took, for a log: pieces, GPU time, rate, the longest piece, the revision
void md_gpu_pieces_describe(const md_gpu_pieces_t* p, char* buf, size_t cap);

#ifdef __cplusplus
}
#endif
