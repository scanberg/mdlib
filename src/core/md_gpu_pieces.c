#include <core/md_gpu_pieces.h>

#include <core/md_common.h>

#include <stdio.h>
#include <string.h>

// A piece's growth over the last one's, at most: one measure of a cheap region does not make the next piece
// ten times larger where the work may be dear
#define PIECES_GROWTH 4.0
// Pieces in flight before the first is measured: a measure soon, and little of unknown cost ahead of the display
#define PIECES_UNMEASURED_IN_FLIGHT 2u
// Polls a passed piece's GPU time is asked for again before the queue is taken not to time
#define PIECES_LATE_POLLS 4u

static uint64_t pieces_round(const md_gpu_pieces_t* p, double want) {
    const uint64_t left = p->total - p->issued;
    uint64_t n = want <= 1.0 ? 1 : (want >= (double)UINT64_MAX / 2 ? UINT64_MAX / 2 : (uint64_t)want);
    n = MD_MAX(n, p->min_units);
    n = MD_MIN(n, p->max_units);
    n = (n / p->min_units) * p->min_units;
    n = MD_MAX(n, p->min_units);
    return MD_MIN(n, left);
}

void md_gpu_pieces_init(md_gpu_pieces_t* p, md_gpu_stream_t stream, uint64_t total, uint64_t first_units, uint64_t min_units, uint64_t max_units,
                        const md_gpu_pieces_policy_t* policy) {
    memset(p, 0, sizeof(*p));
    p->stream    = stream;
    p->total     = total;
    p->min_units = min_units ? min_units : 1;
    p->max_units = max_units ? MD_MAX(max_units, p->min_units) : MD_MAX(total, p->min_units);
    p->first_units = first_units ? first_units : MD_MAX(total / 64, p->min_units);
    p->piece_ms  = policy && policy->piece_ms > 0.0 ? policy->piece_ms : MD_GPU_PIECES_PIECE_MS;
    p->ahead_ms  = policy && policy->ahead_ms != 0.0 ? policy->ahead_ms : MD_GPU_PIECES_AHEAD_MS;
    p->ahead_ms_max = policy && policy->ahead_ms_max > 0.0 ? policy->ahead_ms_max : 4.0 * p->ahead_ms;
}

uint64_t md_gpu_pieces_next(md_gpu_pieces_t* p) {
    if (p->issued >= p->total) return 0;
    const uint32_t in_flight = p->tail - p->head;
    // bounded, a marker per piece; everything at once, the newest marker covers the pieces past the ring's room
    if (in_flight >= MD_GPU_PIECES_MARKS && p->ahead_ms >= 0.0) return 0;
    const bool measured = p->rate > 0.0 && !p->untimed;
    if (p->ahead_ms >= 0.0 && in_flight > 0) {
        if (measured) {
            // the GPU time of what is ahead, by the rate
            if ((double)(p->issued - p->done) / p->rate >= p->ahead_ms) return 0;
        } else if (p->untimed) {
            const double n = p->ahead_ms / p->piece_ms;
            if ((double)in_flight >= MD_MAX(n, 1.0)) return 0;
        } else if (in_flight >= PIECES_UNMEASURED_IN_FLIGHT) {
            return 0;
        }
    }
    double want = (double)p->first_units;
    if (measured) {
        want = p->rate * p->piece_ms;
        if (p->last_units) want = MD_MIN(want, PIECES_GROWTH * (double)p->last_units);
    }
    return pieces_round(p, want);
}

bool md_gpu_pieces_end(md_gpu_pieces_t* p, uint64_t units) {
    units = MD_MIN(units, p->total - p->issued);
    p->issued += units;
    p->last_units = units;
    p->num_pieces += 1;
    const md_gpu_sync_t sync = md_gpu_stream_record(p->stream);
    if (!md_gpu_sync_is_valid(sync)) {
        // nothing was submitted (or ever has been on the stream): nothing to wait for
        p->done = p->issued;
        return units == 0;
    }
    if (p->tail - p->head >= MD_GPU_PIECES_MARKS) {
        // full (next() says so; a caller past it): the newest marker covers this piece too
        md_gpu_pieces_mark_t* m = &p->mark[(p->tail - 1) % MD_GPU_PIECES_MARKS];
        m->sync = sync;
        m->end = p->issued;
        m->units += units;
        return true;
    }
    md_gpu_pieces_mark_t* m = &p->mark[p->tail++ % MD_GPU_PIECES_MARKS];
    m->sync = sync;
    m->end = p->issued;
    m->units = units;
    return true;
}

// A passed piece's GPU time into the rate; false while it is not there
static bool pieces_measure(md_gpu_pieces_t* p, const md_gpu_pieces_mark_t* m) {
    double sec = 0.0;
    if (!md_gpu_sync_gpu_seconds(m->sync, &sec)) return false;
    const double ms = sec * 1.0e3;
    if (ms > 0.0 && m->units > 0) {
        const double sample = (double)m->units / ms;
        // down at once (a dear region: the next piece no longer than the target), up gradually
        if (p->rate <= 0.0 || sample < p->rate) p->rate = sample;
        else p->rate = 0.7 * p->rate + 0.3 * sample;
        p->num_measured += 1;
        p->gpu_ms += ms;
        p->measured_units += m->units;
        p->max_piece_ms = MD_MAX(p->max_piece_ms, ms);
    }
    return true;
}

bool md_gpu_pieces_poll(md_gpu_pieces_t* p) {
    // the passed pieces whose time came late
    uint32_t n = 0;
    for (uint32_t i = 0; i < p->num_late; ++i) {
        if (pieces_measure(p, &p->late[i])) continue;
        const uint8_t polls = (uint8_t)(p->late_polls[i] + 1);
        if (polls < PIECES_LATE_POLLS) {
            p->late[n] = p->late[i];
            p->late_polls[n] = polls;
            n += 1;
        } else if (p->num_measured == 0) {
            p->untimed = true;      // never a time: the queue cannot time
        }
    }
    p->num_late = n;

    const bool any = p->head != p->tail;
    while (p->head != p->tail) {
        md_gpu_pieces_mark_t* m = &p->mark[p->head % MD_GPU_PIECES_MARKS];
        if (!md_gpu_sync_is_complete(m->sync)) break;
        p->done = m->end;
        if (!p->untimed && !pieces_measure(p, m)) {
            if (p->num_late < MD_COUNTOF(p->late)) {
                p->late[p->num_late] = *m;
                p->late_polls[p->num_late] = 0;
                p->num_late += 1;
            }
        }
        p->head += 1;
    }
    // every piece passed with more to issue: the GPU had nothing to do since some time before this poll
    if (any && p->head == p->tail && p->issued < p->total && p->ahead_ms >= 0.0 && (p->rate > 0.0 || p->untimed)) {
        p->num_dry += 1;
        p->ahead_ms = MD_MIN(2.0 * p->ahead_ms, p->ahead_ms_max);
    }
    return p->done >= p->total;
}

void md_gpu_pieces_wait(md_gpu_pieces_t* p) {
    if (p->head != p->tail) md_gpu_sync_wait(p->mark[p->head % MD_GPU_PIECES_MARKS].sync);
}

void md_gpu_pieces_describe(const md_gpu_pieces_t* p, char* buf, size_t cap) {
    if (!buf || !cap) return;
    snprintf(buf, cap, "%u pieces, %.1f ms of GPU time over %u measured (%.3g units/ms, the longest %.1f ms), %.0f ms ahead, %u dry polls%s; " MD_GPU_PIECES_REVISION,
             p->num_pieces, p->gpu_ms, p->num_measured, p->measured_units && p->gpu_ms > 0.0 ? (double)p->measured_units / p->gpu_ms : 0.0, p->max_piece_ms,
             p->ahead_ms, p->num_dry, p->untimed ? ", untimed" : "");
}
