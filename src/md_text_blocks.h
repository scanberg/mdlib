#pragma once

// For the loaders of text formats which read their lines in parallel (md_executor.h): a file, or a string, a block at a
// time, and a block in pieces at line boundaries. Internal to mdlib.
//
// A block is whole lines: what was read, up to its last newline. The rest is kept for the next block, so a line is
// never split between two. A piece is a range of a block which begins at the start of a line and ends after a newline
// (or at the end of the block).

#include <core/md_common.h>
#include <core/md_str.h>
#include <core/md_os.h>
#include <core/md_allocator.h>

#include <string.h>

#define MD_TEXT_BLOCK_SIZE      MEGABYTES(32)
#define MD_TEXT_MIN_PIECE_SIZE  KILOBYTES(256)

typedef struct md_text_source_t {
    md_file_t       file;       // when is_file; a string otherwise, and then data is all of it
    bool            is_file;
    char*           buf;
    size_t          cap;
    const char*     data;
    size_t          len;
    size_t          pos;        // what is read of data
    bool            eof;
    size_t          total;      // bytes of all of it: the file's size, the string's length
    md_allocator_i* alloc;
} md_text_source_t;

// Keeps what is unread, and reads as much as fits after it
static inline void md_text_source_fill(md_text_source_t* src) {
    if (!src->is_file || src->eof) return;
    const size_t rest = src->len - src->pos;
    if (rest == src->cap) return;   // Full: nothing read is taken yet
    if (rest && src->pos) memmove(src->buf, src->buf + src->pos, rest);
    src->pos = 0;
    src->len = rest;
    const size_t bytes = md_file_read(src->file, src->buf + rest, src->cap - rest);
    src->len += bytes;
    // Short of what was asked for: the end of the file
    if (src->len < src->cap) src->eof = true;
}

static inline md_text_source_t md_text_source_from_str(str_t str) {
    md_text_source_t src = { .data = str.ptr, .len = str.len, .eof = true, .total = str.len };
    return src;
}

static inline bool md_text_source_open(md_text_source_t* src, str_t filename, md_allocator_i* alloc) {
    MEMSET(src, 0, sizeof(*src));
    if (!md_file_open(&src->file, filename, MD_FILE_READ)) return false;
    src->is_file = true;
    src->alloc = alloc;
    src->total = (size_t)md_file_size(src->file);
    src->cap = MD_TEXT_BLOCK_SIZE;
    src->buf = (char*)md_alloc(alloc, src->cap);
    src->data = src->buf;
    md_text_source_fill(src);
    return true;
}

static inline void md_text_source_close(md_text_source_t* src) {
    if (src->is_file) {
        md_file_close(&src->file);
        md_free(src->alloc, src->buf, src->cap);
    }
    MEMSET(src, 0, sizeof(*src));
}

// The next line, without its newline (nor a carriage return before it). False at the end.
static inline bool md_text_source_line(str_t* line, md_text_source_t* src) {
    for (;;) {
        const char* beg = src->data + src->pos;
        const char* nl  = (const char*)memchr(beg, '\n', src->len - src->pos);
        if (nl) {
            *line = (str_t){ beg, (size_t)(nl - beg) };
            src->pos += (size_t)(nl - beg) + 1;
            break;
        }
        if (src->eof || (src->pos == 0 && src->len == src->cap)) {
            // The last line, without a newline; or one longer than a block, cut at it
            if (src->pos == src->len) return false;
            *line = (str_t){ beg, src->len - src->pos };
            src->pos = src->len;
            break;
        }
        md_text_source_fill(src);
    }
    if (line->len && line->ptr[line->len - 1] == '\r') line->len -= 1;
    return true;
}

// The next line without taking it
static inline str_t md_text_source_peek(md_text_source_t* src) {
    if (src->pos == src->len) md_text_source_fill(src);
    const char* beg = src->data + src->pos;
    const char* nl  = (const char*)memchr(beg, '\n', src->len - src->pos);
    str_t line = { beg, nl ? (size_t)(nl - beg) : src->len - src->pos };
    if (line.len && line.ptr[line.len - 1] == '\r') line.len -= 1;
    return line;
}

// The next block of whole lines, read as far as it goes: empty at the end. Taken with md_text_source_advance.
static inline str_t md_text_source_block(md_text_source_t* src) {
    if (src->is_file && !src->eof) md_text_source_fill(src);
    const char* text = src->data + src->pos;
    size_t len = src->len - src->pos;
    if (!src->eof) {
        const char* c = text + len;
        while (c > text && c[-1] != '\n') --c;
        // A line longer than a block: cut at it
        if (c > text) len = (size_t)(c - text);
    }
    return (str_t){ text, len };
}

static inline void md_text_source_advance(md_text_source_t* src, size_t bytes) {
    src->pos += bytes;
}

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
    str_t line = { c, (size_t)(end - c) };
    if (line.len && line.ptr[line.len - 1] == '\r') line.len -= 1;
    return line;
}

// The line at c, up to e, and where the next begins
static inline str_t md_text_next_line(const char** c, const char* e) {
    const char* beg = *c;
    const char* nl = (const char*)memchr(beg, '\n', (size_t)(e - beg));
    const char* le = nl ? nl : e;
    *c = nl ? nl + 1 : e;
    str_t line = { beg, (size_t)(le - beg) };
    if (line.len && line.ptr[line.len - 1] == '\r') line.len -= 1;
    return line;
}
