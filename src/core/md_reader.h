#pragma once

#include <core/md_common.h>
#include <core/md_str.h>
#include <core/md_os.h>     // md_file_t
#include <core/md_intrinsics.h>

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

// [beg, end[ without a carriage return at its end
static inline str_t md_reader_strip_cr_raw(const char* beg, const char* end) {
    if (end > beg && end[-1] == '\r') --end;
    return (str_t){ beg, (size_t)(end - beg) };
}

// ### FINDING THE ENDS OF LINES ###
//
// md_text_find_newline: the first '\n' of [c, e), else e. What every line here is found by, so in place of memchr
// (a call, and a setup for any length) a search inline, a register of bytes at a time: compared with '\n' at once and
// the matches as a mask, one bit per byte (SSE2, AVX2), four per byte on arm64, where the narrowing shift of the
// comparison (shrn) is what movemask is on x86.
//
// The loads are aligned, so none crosses a page: past e (or before c, in the register of c) is read only where a page
// holds bytes of [c, e) as well, which may be read. Address sanitizers see the reads past the range: they are told not
// to look.

#if defined(__AVX2__)
#   include <immintrin.h>
#   define MD_TEXT_SCAN_WIDTH 32
#   define MD_TEXT_SCAN_BITS  1
#elif defined(__x86_64__) || defined(_M_X64) || defined(__SSE2__)
#   include <emmintrin.h>
#   define MD_TEXT_SCAN_WIDTH 16
#   define MD_TEXT_SCAN_BITS  1
#elif defined(__aarch64__) || defined(_M_ARM64)
#   include <arm_neon.h>
#   define MD_TEXT_SCAN_WIDTH 16
#   define MD_TEXT_SCAN_BITS  4
#endif

#if defined(__clang__) || defined(__GNUC__)
#   define MD_TEXT_NO_SANITIZE __attribute__((no_sanitize_address))
#elif defined(_MSC_VER)
#   define MD_TEXT_NO_SANITIZE __declspec(no_sanitize_address)
#else
#   define MD_TEXT_NO_SANITIZE
#endif

#if defined(MD_TEXT_SCAN_WIDTH)
// The newlines of the aligned register of bytes at p: MD_TEXT_SCAN_BITS bits of the mask per byte, the first byte lowest
MD_TEXT_NO_SANITIZE static inline uint64_t md_text_newline_mask(const char* p) {
#   if defined(__AVX2__)
    return (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_load_si256((const __m256i*)p), _mm256_set1_epi8('\n')));
#   elif MD_TEXT_SCAN_BITS == 1
    return (uint32_t)_mm_movemask_epi8(_mm_cmpeq_epi8(_mm_load_si128((const __m128i*)p), _mm_set1_epi8('\n')));
#   else
    const uint8x16_t eq = vceqq_u8(vld1q_u8((const uint8_t*)p), vdupq_n_u8('\n'));
    return vget_lane_u64(vreinterpret_u64_u8(vshrn_n_u16(vreinterpretq_u16_u8(eq), 4)), 0);
#   endif
}
#endif

// THE LINES OF A RANGE, in bulk: md_text_lines_t takes the lines of [c, e) one after the other, from masks of the
// newlines of a chunk of bytes at a time (64 bytes, 16 on arm64), each line from the mask alone (its lowest bit,
// cleared). Where md_text_find_newline searches from where the last line ended, a chain of loads which each wait for
// the one before, the chunks are loaded one after the other whatever the lines are: twice as fast over lines of 80
// characters (8 GB/s with AVX2, 6 with SSE2, against 4.4 for memchr or md_text_find_newline line by line).

#if defined(MD_TEXT_SCAN_WIDTH)
#   define MD_TEXT_CHUNK (64 / MD_TEXT_SCAN_BITS)   // Bytes of a chunk: one bit of a 64 bit mask per byte, four on arm64

// The newlines of the aligned chunk at p, MD_TEXT_SCAN_BITS bits per byte
MD_TEXT_NO_SANITIZE static inline uint64_t md_text_newline_chunk(const char* p) {
#   if defined(__AVX2__)
    return md_text_newline_mask(p) | (md_text_newline_mask(p + 32) << 32);
#   elif MD_TEXT_SCAN_BITS == 1
    return md_text_newline_mask(p) | (md_text_newline_mask(p + 16) << 16) | (md_text_newline_mask(p + 32) << 32) | (md_text_newline_mask(p + 48) << 48);
#   else
    return md_text_newline_mask(p);
#   endif
}
#endif

typedef struct md_text_lines_t {
    const char* beg;        // Of the next line
    const char* end;        // Of the range
#if defined(MD_TEXT_SCAN_WIDTH)
    const char* chunk;      // The aligned chunk of mask
    uint64_t    mask;       // Its newlines not taken yet
#endif
} md_text_lines_t;

MD_TEXT_NO_SANITIZE static inline md_text_lines_t md_text_lines(const char* c, const char* e) {
    md_text_lines_t it = { .beg = c, .end = e };
#if defined(MD_TEXT_SCAN_WIDTH)
    it.chunk = (const char*)((uintptr_t)c & ~(uintptr_t)(MD_TEXT_CHUNK - 1));
    it.mask = c < e ? md_text_newline_chunk(it.chunk) & (~UINT64_C(0) << ((uintptr_t)(c - it.chunk) * MD_TEXT_SCAN_BITS)) : 0;
#endif
    return it;
}

// The next line, without its newline (nor a carriage return before it): false after the last. The last line of the
// range need not end in a newline.
MD_TEXT_NO_SANITIZE static inline bool md_text_lines_next(md_text_lines_t* it, str_t* line) {
    const char* nl = it->end;
#if defined(MD_TEXT_SCAN_WIDTH)
    for (;;) {
        if (it->mask) {
            const uint32_t bit = (uint32_t)ctz64(it->mask);
            // The lowest newline: its bit, or its nibble on arm64 (whose bits are all set)
            it->mask &= ~(((UINT64_C(1) << MD_TEXT_SCAN_BITS) - 1) << bit);
            nl = it->chunk + bit / MD_TEXT_SCAN_BITS;
            if (nl > it->end) nl = it->end;
            break;
        }
        it->chunk += MD_TEXT_CHUNK;
        if (it->chunk >= it->end) break;
        it->mask = md_text_newline_chunk(it->chunk);
    }
#else
    const char* found = (const char*)memchr(it->beg, '\n', (size_t)(it->end - it->beg));
    if (found) nl = found;
#endif
    if (it->beg >= it->end) return false;
    *line = md_reader_strip_cr_raw(it->beg, nl);
    it->beg = nl < it->end ? nl + 1 : it->end;
#if defined(MD_TEXT_SCAN_WIDTH)
    if (it->beg >= it->end) it->mask = 0;
#endif
    return true;
}

// The newlines of [c, e): a count of the bits of the masks of its chunks
MD_TEXT_NO_SANITIZE static inline size_t md_text_count_newlines(const char* c, const char* e) {
    if (c >= e) return 0;
#if defined(MD_TEXT_SCAN_WIDTH)
    const char* chunk = (const char*)((uintptr_t)c & ~(uintptr_t)(MD_TEXT_CHUNK - 1));
    uint64_t m = md_text_newline_chunk(chunk) & (~UINT64_C(0) << ((uintptr_t)(c - chunk) * MD_TEXT_SCAN_BITS));
    size_t n = 0;
    for (;;) {
        const size_t past = (size_t)MAX((ptrdiff_t)0, (chunk + MD_TEXT_CHUNK) - e);
        if (past) {
            // The last chunk: not past e
            m &= ~UINT64_C(0) >> (past * MD_TEXT_SCAN_BITS);
            return n + (size_t)popcnt64(m) / MD_TEXT_SCAN_BITS;
        }
        n += (size_t)popcnt64(m) / MD_TEXT_SCAN_BITS;
        chunk += MD_TEXT_CHUNK;
        if (chunk >= e) return n;   // e at the end of a chunk: none past it is loaded
        m = md_text_newline_chunk(chunk);
    }
#else
    size_t n = 0;
    while (c < e && (c = (const char*)memchr(c, '\n', (size_t)(e - c))) != NULL) { ++n; ++c; }
    return n;
#endif
}

MD_TEXT_NO_SANITIZE static inline const char* md_text_find_newline(const char* c, const char* e) {
#if defined(MD_TEXT_SCAN_WIDTH)
    if (c >= e) return e;
    const char* a = (const char*)((uintptr_t)c & ~(uintptr_t)(MD_TEXT_SCAN_WIDTH - 1));
    uint64_t m = md_text_newline_mask(a) >> ((uintptr_t)(c - a) * MD_TEXT_SCAN_BITS);
    if (m) {
        const char* nl = c + ctz64(m) / MD_TEXT_SCAN_BITS;
        return nl < e ? nl : e;
    }
    for (a += MD_TEXT_SCAN_WIDTH; a < e; a += MD_TEXT_SCAN_WIDTH) {
        m = md_text_newline_mask(a);
        if (m) {
            const char* nl = a + ctz64(m) / MD_TEXT_SCAN_BITS;
            return nl < e ? nl : e;
        }
    }
    return e;
#else
    const char* nl = (const char*)memchr(c, '\n', (size_t)(e - c));
    return nl ? nl : e;
#endif
}

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
    const char* nl = md_text_find_newline(r->cur, r->end);
    if (nl < r->end) {
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
#define MD_TEXT_PASS_SIZE       MEGABYTES(4)    // What a loader reads in two passes (lines, then fields): in cache for the second

// The whole lines of the first MD_TEXT_PASS_SIZE bytes of a block of whole lines, all of it when it is no longer (or the
// first line is longer): the loader advances by what it took, and the next block is the rest
static inline size_t md_text_pass_len(const char* text, size_t len) {
    if (len <= MD_TEXT_PASS_SIZE) return len;
    size_t at = MD_TEXT_PASS_SIZE;
    while (at > 0 && text[at - 1] != '\n') --at;
    return at > 0 ? at : len;
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
            const char* nl = md_text_find_newline(text + e, text + len);
            e = nl < text + len ? (size_t)(nl - text) + 1 : len;
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
    const char* nl = md_text_find_newline(beg, e);
    const char* le = nl;
    *c = nl < e ? nl + 1 : e;
    return md_reader_strip_cr((str_t){ beg, (size_t)(le - beg) });
}

#ifdef __cplusplus
}
#endif
