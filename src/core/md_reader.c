#include <core/md_reader.h>

#include <core/md_allocator.h>
#include <core/md_log.h>

#include <stddef.h>
#include <string.h>

// The readers are handed over as &x.r, and their vtables cast r back to the reader
STATIC_ASSERT(offsetof(md_file_reader_t, r) == 0, "md_reader_i must be the first member");

// ### THE INTERFACE ###

bool md_reader_init(md_reader_i* r, const md_reader_vtable_t* vt, size_t cap, int64_t size, struct md_allocator_i* alloc) {
    ASSERT(r);
    ASSERT(vt && vt->read);
    ASSERT(alloc);
    MEMSET(r, 0, sizeof(*r));
    if (cap == 0) cap = MD_READER_DEFAULT_CAPACITY;
    r->mem = (char*)md_alloc(alloc, cap + MD_READER_PADDING);
    if (!r->mem) {
        MD_LOG_ERROR("md_reader: out of memory for a buffer of %zu bytes", cap);
        return false;
    }
    // The padding is read, never written: what it holds is never part of a line, but it is defined
    MEMSET(r->mem + cap, 0, MD_READER_PADDING);
    r->buf = r->cur = r->end = r->mem;
    r->cap = cap;
    r->vt = vt;
    r->size = size;
    r->alloc = alloc;
    return true;
}

void md_reader_free(md_reader_i* r) {
    ASSERT(r);
    if (r->mem) md_free(r->alloc, r->mem, r->cap + MD_READER_PADDING);
    MEMSET(r, 0, sizeof(*r));
}

// Keeps what is not taken at the start of the buffer, and returns the room after it: 0 at the end of the source, for a
// string, or when the buffer is full of what is not taken
static size_t refill_begin(md_reader_i* r) {
    if (!r->vt || (r->flags & MD_READER_FLAG_END)) return 0;
    const size_t rest = (size_t)(r->end - r->cur);
    if (rest == r->cap) return 0;
    r->base += (int64_t)(r->cur - r->buf);
    if (rest && r->cur != r->mem) memmove(r->mem, r->cur, rest);
    r->buf = r->cur = r->mem;
    r->end = r->mem + rest;
    return r->cap - rest;
}

// n of the want bytes of room were read
static void refill_end(md_reader_i* r, size_t n, size_t want) {
    r->end += n;
    if (n < want) r->flags |= MD_READER_FLAG_END;
}

// Keeps what is not taken at the start of the buffer, and reads as much after it as fits. Returns the bytes read: 0 at
// the end of the source, for a string, or when the buffer is full of what is not taken.
static size_t refill(md_reader_i* r) {
    const size_t want = refill_begin(r);
    if (want == 0) return 0;
    const size_t n = r->vt->read(r, r->mem + (r->end - r->mem), want);
    refill_end(r, n, want);
    return n;
}

bool md_reader_line_slow(md_reader_i* r, str_t* line, const char** next) {
    // [cur, cur + searched[ has no newline: what was there before reading more
    size_t searched = (size_t)(r->end - r->cur);
    for (;;) {
        if (!refill(r)) {
            // The last line, without a newline; or one longer than the buffer, cut at it
            if (r->cur == r->end) return false;
            *line = md_reader_strip_cr((str_t){ r->cur, (size_t)(r->end - r->cur) });
            *next = r->end;
            return true;
        }
        const char* from = r->cur + searched;
        const char* nl = md_text_find_newline(from, r->end);
        if (nl < r->end) {
            *line = md_reader_strip_cr((str_t){ r->cur, (size_t)(nl - r->cur) });
            *next = nl + 1;
            return true;
        }
        searched = (size_t)(r->end - r->cur);
    }
}

// The whole lines of what the buffer holds: up to the last newline (the rest is read again with the next block)
static str_t whole_lines(const md_reader_i* r) {
    const char* text = r->cur;
    size_t len = (size_t)(r->end - r->cur);
    if (!(r->flags & MD_READER_FLAG_END)) {
        const char* c = text + len;
        while (c > text && c[-1] != '\n') --c;
        if (c > text) len = (size_t)(c - text);     // else a line longer than the buffer: cut at it
    }
    return (str_t){ text, len };
}

str_t md_reader_block(md_reader_i* r) {
    ASSERT(r);
    refill(r);
    return whole_lines(r);
}

// The reads of a refill over the threads of an executor: pieces of the room, each read at its offset
#define READER_PARALLEL_PIECE MEGABYTES(2)

typedef struct parallel_read_t {
    const md_reader_i* r;
    int64_t offset;     // Of the source, where the room begins
    char*   dst;
    size_t  want;
    size_t* got;        // By piece
} parallel_read_t;

static void parallel_read_piece(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    parallel_read_t* pr = (parallel_read_t*)arg;
    for (uint32_t p = beg; p < end; ++p) {
        const size_t at = (size_t)p * READER_PARALLEL_PIECE;
        const size_t len = MIN(READER_PARALLEL_PIECE, pr->want - at);
        pr->got[p] = pr->r->vt->read_at(pr->r, pr->offset + (int64_t)at, pr->dst + at, len);
    }
}

str_t md_reader_block_parallel(md_reader_i* r, const md_executor_t* ex) {
    ASSERT(r);
    if (md_executor_num_threads(ex) <= 1 || !r->vt || !r->vt->read_at || !r->vt->seek) {
        return md_reader_block(r);
    }
    const size_t want = refill_begin(r);
    if (want > 0) {
        const uint32_t num_pieces = (uint32_t)((want + READER_PARALLEL_PIECE - 1) / READER_PARALLEL_PIECE);
        size_t got_buf[64];
        size_t* got = num_pieces <= 64 ? got_buf : (size_t*)md_alloc(r->alloc, num_pieces * sizeof(size_t));
        parallel_read_t pr = {
            .r = r,
            .offset = r->base + (int64_t)(r->end - r->buf),
            .dst = r->mem + (r->end - r->mem),
            .want = want,
            .got = got,
        };
        // Cancelled: nothing more is read (the pieces not read count as none)
        for (uint32_t p = 0; p < num_pieces; ++p) got[p] = 0;
        md_parallel_for(ex, num_pieces, 1, parallel_read_piece, &pr);
        // What was read without a gap: a piece short of its length is the end of the source
        size_t n = 0;
        for (uint32_t p = 0; p < num_pieces; ++p) {
            const size_t len = MIN(READER_PARALLEL_PIECE, want - (size_t)p * READER_PARALLEL_PIECE);
            n += got[p];
            if (got[p] < len) break;
        }
        if (got != got_buf) md_free(r->alloc, got, num_pieces * sizeof(size_t));
        // Where the next read begins, past what was read here
        if (!r->vt->seek(r, pr.offset + (int64_t)n)) r->flags |= MD_READER_FLAG_ERROR;
        refill_end(r, n, want);
    }
    return whole_lines(r);
}

bool md_reader_seek(md_reader_i* r, int64_t offset) {
    ASSERT(r);
    if (offset < 0) return false;
    // Within what is buffered (all of a string)
    if (offset >= r->base && offset <= r->base + (int64_t)(r->end - r->buf)) {
        r->cur = r->buf + (offset - r->base);
        return true;
    }
    if (!r->vt || !r->vt->seek) return false;
    if (r->size >= 0 && offset > r->size) return false;
    if (!r->vt->seek(r, offset)) return false;
    r->base = offset;
    r->buf = r->cur = r->end = r->mem;
    r->flags &= ~(uint32_t)MD_READER_FLAG_END;
    return true;
}

// ### FILE READER ###

// A file opened with MD_FILE_CONCURRENT (md_file_reader_open) has no position: the reader keeps its own (pos) and reads
// at it. One opened by the caller (md_file_reader_init) is read from its position, which reading moves.
static bool file_concurrent(const md_file_reader_t* fr) {
    return (fr->file.flags & MD_FILE_CONCURRENT) != 0;
}

static size_t file_read(md_reader_i* r, char* dst, size_t cap) {
    md_file_reader_t* fr = (md_file_reader_t*)r;
    const bool concurrent = file_concurrent(fr);
    size_t total = 0;
    // A read may come back short of what was asked for before the end (large requests on some platforms)
    while (total < cap) {
        const size_t n = concurrent ? md_file_read_at(fr->file, fr->pos, dst + total, cap - total)
                                    : md_file_read(fr->file, dst + total, cap - total);
        if (n == 0) break;
        total += n;
        fr->pos += (int64_t)n;
    }
    return total;
}

static bool file_seek(md_reader_i* r, int64_t offset) {
    md_file_reader_t* fr = (md_file_reader_t*)r;
    if (!file_concurrent(fr) && !md_file_seek(fr->file, offset, MD_FILE_BEG)) return false;
    fr->pos = offset;
    return true;
}

static size_t file_read_at(const md_reader_i* r, int64_t offset, char* dst, size_t cap) {
    const md_file_reader_t* fr = (const md_file_reader_t*)r;
    return md_file_read_at(fr->file, offset, dst, cap);
}

static const md_reader_vtable_t file_vtable = {
    .read = file_read,
    .seek = file_seek,
    .read_at = file_read_at,
    .name = "file",
};

bool md_file_reader_init(md_file_reader_t* fr, md_file_t file, size_t cap, struct md_allocator_i* alloc) {
    ASSERT(fr);
    if (!md_file_valid(file)) {
        MEMSET(fr, 0, sizeof(*fr));
        return false;
    }
    if (!md_reader_init(&fr->r, &file_vtable, cap, (int64_t)md_file_size(file), alloc)) {
        return false;
    }
    fr->file = file;
    fr->owns_file = false;
    // From where the file is positioned: the start of one which has no position
    fr->pos = (file.flags & MD_FILE_CONCURRENT) ? 0 : (int64_t)md_file_tell(file);
    fr->r.base = fr->pos;
    return true;
}

bool md_file_reader_open(md_file_reader_t* fr, str_t path, size_t cap, struct md_allocator_i* alloc) {
    ASSERT(fr);
    // Read front to back, and at offsets from several threads (md_reader_block_parallel): a handle every thread reads
    // at once, which on Windows a synchronous one would serialize
    md_file_t file = {0};
    if (!md_file_open(&file, path, MD_FILE_READ | MD_FILE_CONCURRENT | MD_FILE_SEQUENTIAL)) {
        MEMSET(fr, 0, sizeof(*fr));
        return false;
    }
    if (!md_file_reader_init(fr, file, cap, alloc)) {
        md_file_close(&file);
        return false;
    }
    fr->owns_file = true;
    return true;
}

void md_file_reader_close(md_file_reader_t* fr) {
    ASSERT(fr);
    if (fr->owns_file) md_file_close(&fr->file);
    md_reader_free(&fr->r);
    MEMSET(fr, 0, sizeof(*fr));
}
