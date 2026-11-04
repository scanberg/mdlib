#pragma once

#include <core/md_common.h>
#include <core/md_str.h>
#include <core/md_os.h>     // md_file_t

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

struct md_allocator_i;

/*

READERS

Where input comes from: the counterpart of the writers (md_writer.h). A parser reads its text through a reader
interface (md_reader_i), and its caller decides where the text comes from by which reader it hands over:

    md_reader_from_str()    a string in memory, all of it at once
    md_file_reader_t        a file, through a buffer of its own

    md_file_reader_t fr;
    if (md_file_reader_open(&fr, filename, MEGABYTES(1), alloc)) {
        str_t line;
        while (md_reader_extract_line(&fr.r, &line)) { ... }
        md_file_reader_close(&fr);
    }

    md_reader_i sr = md_reader_from_str(text);
    parse(&sr);

THE INTERFACE

Every reader has an md_reader_i as its first member, r, and &x.r is what is handed over. The interface holds what is
read and not yet taken, [cur, end[: taking a line is a search in it, inline. Only a read which runs off the end of it
consults the kind of reader, through its vtable, to read more. A string reader has no vtable: all of it is there from
the start. A source of another kind (a compressed file) is a vtable of its own.

A line is a view into the reader's buffer (or the string): valid until the next call which reads more, which is any
call taking or peeking a line or a block. One copies what one keeps.

A line is its text without its newline, and without a carriage return before it (CRLF files). The last line of a source
need not end in a newline. A line longer than the buffer is cut at its length (and its remainder is the next line):
buffers are made larger than any line of the formats read through them.

BLOCKS

A loader which parses its lines in parallel (md_executor.h) takes them a block at a time: md_reader_block gives the
whole lines the buffer holds (reading more first), md_reader_advance takes what was parsed. md_text_split cuts a block
into pieces at line boundaries, one per task.

PADDING

The buffer of a file reader has MD_READER_PADDING bytes past its end, which may be read (not written): a search over
a block may load whole registers without checking for its end. A string reader has none, as the string is the caller's.

POSITION

md_reader_tell is the offset in the source of the next byte to take, and md_reader_seek goes to one: anywhere in a
string or a file, within what is buffered for a source which cannot seek (a stream). An offset taken with
md_reader_tell before a line is taken is the line's: what the trajectory loaders record to read a frame again.

A reader is used by one thread at a time.

*/

#define MD_READER_PADDING 64

// The default size of a file reader's buffer: lines of any format we read fit many times over
#define MD_READER_DEFAULT_CAPACITY MEGABYTES(1)

typedef struct md_reader_i md_reader_i;

typedef enum md_reader_flags_t {
    MD_READER_FLAG_END   = 1u << 0,     // The source is read to its end: what is left is [cur, end[
    MD_READER_FLAG_ERROR = 1u << 1,     // The source failed (a read error). It sticks, and ends the reading
} md_reader_flags_t;

// What makes a kind of reader. Implemented for each source, and for a source of one's own.
typedef struct md_reader_vtable_t {
    // Reads up to cap bytes into dst and returns how many: fewer only at the end of the source (or on failure, with
    // MD_READER_FLAG_ERROR set).
    size_t (*read)(md_reader_i* r, char* dst, size_t cap);

    // Positions the source at offset from its start, where the next read begins. NULL: it cannot (a stream).
    bool (*seek)(md_reader_i* r, int64_t offset);

    const char* name;   // "file", ...: for the debugger
} md_reader_vtable_t;

struct md_reader_i {
    const char* buf;    // What is read, [buf, end[: [buf, cur[ is taken, [cur, end[ not yet
    const char* cur;
    const char* end;

    const md_reader_vtable_t* vt;   // NULL: a string, all of which is in [buf, end[
    char*    mem;                   // The buffer a reader of a source reads into: cap bytes, and the padding after
    size_t   cap;
    int64_t  base;                  // The offset in the source of buf[0]
    int64_t  size;                  // The bytes of the whole source, -1 when it is not known (a stream)
    struct md_allocator_i* alloc;   // What mem is allocated with
    uint32_t flags;                 // md_reader_flags_t
};

// ### CONCRETE READERS ###

typedef struct md_file_reader_t {
    md_reader_i r;
    md_file_t file;
    bool owns_file;     // Opened by md_file_reader_open, closed by md_file_reader_close
} md_file_reader_t;

#ifdef __cplusplus
extern "C" {
#endif

// ### THE INTERFACE ###

// A reader of a string: all of it, at once. Nothing to free.
static inline md_reader_i md_reader_from_str(str_t str) {
    md_reader_i r = {0};
    r.buf = r.cur = str.ptr;
    r.end = str.ptr + str.len;
    r.size = (int64_t)str.len;
    r.flags = MD_READER_FLAG_END;
    return r;
}

// For a reader of one's own: a struct with an md_reader_i as its first member, whose vtable casts r back to it. Its
// buffer, of cap bytes (and the padding), is allocated with alloc. size is the bytes of the source, -1 if not known.
bool md_reader_init(md_reader_i* r, const md_reader_vtable_t* vt, size_t cap, int64_t size, struct md_allocator_i* alloc);

// Frees the buffer of a reader made with md_reader_init (a string reader has none)
void md_reader_free(md_reader_i* r);

// The parts of the line functions which run off the end of what is read. Not to be called directly.
bool md_reader_line_slow(md_reader_i* r, str_t* line, const char** next);

static inline str_t md_reader_strip_cr(str_t line) {
    if (line.len && line.ptr[line.len - 1] == '\r') line.len -= 1;
    return line;
}

// The next line, in line, and where the one after it begins: false at the end
static inline bool md_reader_find_line(md_reader_i* r, str_t* line, const char** next) {
    const char* nl = (const char*)memchr(r->cur, '\n', (size_t)(r->end - r->cur));
    if (nl) {
        *line = md_reader_strip_cr((str_t){ r->cur, (size_t)(nl - r->cur) });
        *next = nl + 1;
        return true;
    }
    return md_reader_line_slow(r, line, next);
}

// The next line, taken. False at the end of the source.
static inline bool md_reader_extract_line(md_reader_i* r, str_t* line) {
    ASSERT(r && line);
    const char* next;
    if (!md_reader_find_line(r, line, &next)) return false;
    r->cur = next;
    return true;
}

// The next line, not taken
static inline bool md_reader_peek_line(md_reader_i* r, str_t* line) {
    ASSERT(r && line);
    const char* next;
    return md_reader_find_line(r, line, &next);
}

static inline bool md_reader_skip_line(md_reader_i* r) {
    ASSERT(r);
    str_t line;
    const char* next;
    if (!md_reader_find_line(r, &line, &next)) return false;
    r->cur = next;
    return true;
}

// The whole lines which the buffer holds, as far as it reads, after reading more: empty at the end. A line longer than
// the buffer is cut at it. Taken with md_reader_advance, up to its length.
str_t md_reader_block(md_reader_i* r);

static inline void md_reader_advance(md_reader_i* r, size_t bytes) {
    ASSERT(r && bytes <= (size_t)(r->end - r->cur));
    r->cur += bytes;
}

// The offset in the source of the next byte to take
static inline int64_t md_reader_tell(const md_reader_i* r) {
    return r->base + (int64_t)(r->cur - r->buf);
}

// To the byte of the source at offset. False if the reader cannot go there (a stream, which cannot seek, can still go
// back within what it has buffered; an offset past the end).
bool md_reader_seek(md_reader_i* r, int64_t offset);

static inline bool md_reader_rewind(md_reader_i* r) {
    return md_reader_seek(r, 0);
}

// The bytes of the whole source, -1 when it is not known
static inline int64_t md_reader_size(const md_reader_i* r) {
    return r->size;
}

// Everything is taken
static inline bool md_reader_end(const md_reader_i* r) {
    return r->cur == r->end && (r->flags & MD_READER_FLAG_END);
}

static inline bool md_reader_failed(const md_reader_i* r) {
    return (r->flags & MD_READER_FLAG_ERROR) != 0;
}

// ### FILE READER ###

// Of the file at path, opened for reading, through a buffer of cap bytes (0: MD_READER_DEFAULT_CAPACITY) allocated with
// alloc. False if the file cannot be opened.
bool md_file_reader_open(md_file_reader_t* fr, str_t path, size_t cap, struct md_allocator_i* alloc);

// Of a file opened by the caller, from where it is positioned, which the caller closes after md_file_reader_close
bool md_file_reader_init(md_file_reader_t* fr, md_file_t file, size_t cap, struct md_allocator_i* alloc);

// Frees the buffer, and closes the file if the reader opened it
void md_file_reader_close(md_file_reader_t* fr);

// ### BLOCKS AND PIECES ###
// For the loaders which parse the lines of a block in parallel: a block cut into pieces at line boundaries, one per
// task, of which the lines are taken with md_text_next_line.

#define MD_TEXT_BLOCK_SIZE      MEGABYTES(32)   // The buffer of a loader which parses a block at a time
#define MD_TEXT_MIN_PIECE_SIZE  KILOBYTES(256)

typedef struct md_text_piece_t {
    size_t beg, end;    // bytes of the block
} md_text_piece_t;

// How many pieces to cut a block of len bytes into, for an executor of num_threads: a few per thread to even out the
// lines that cost more, none smaller than MD_TEXT_MIN_PIECE_SIZE (but the one of a small block)
static inline size_t md_text_num_pieces(size_t len, uint32_t num_threads) {
    if (num_threads <= 1) return 1;
    const size_t by_size = len / MD_TEXT_MIN_PIECE_SIZE;
    return MAX((size_t)1, MIN((size_t)num_threads * 4, by_size));
}

// The pieces of a block of whole lines, at line boundaries: some may be empty when lines are long
static inline void md_text_split(md_text_piece_t* pieces, size_t num_pieces, const char* text, size_t len) {
    size_t at = 0;
    for (size_t p = 0; p < num_pieces; ++p) {
        size_t e = (p + 1 == num_pieces) ? len : MAX(at, len * (p + 1) / num_pieces);
        if (e < len && e > 0 && text[e - 1] != '\n') {
            const char* nl = (const char*)memchr(text + e, '\n', len - e);
            e = nl ? (size_t)(nl - text) + 1 : len;
        }
        pieces[p] = (md_text_piece_t){ at, e };
        at = e;
    }
}

// The line which ends right before offset at (a line start), without its newline; at must be > 0
static inline str_t md_text_line_before(const char* text, size_t at) {
    const char* end = text + at - 1;    // its newline
    const char* c = end;
    while (c > text && c[-1] != '\n') --c;
    return md_reader_strip_cr((str_t){ c, (size_t)(end - c) });
}

// The line at c, up to e, and where the next begins
static inline str_t md_text_next_line(const char** c, const char* e) {
    const char* beg = *c;
    const char* nl = (const char*)memchr(beg, '\n', (size_t)(e - beg));
    const char* le = nl ? nl : e;
    *c = nl ? nl + 1 : e;
    return md_reader_strip_cr((str_t){ beg, (size_t)(le - beg) });
}

#ifdef __cplusplus
}
#endif
