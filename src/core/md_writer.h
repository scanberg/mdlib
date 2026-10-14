#pragma once

#include <core/md_common.h>
#include <core/md_str.h>
#include <core/md_os.h>     // md_file_t
#include <core/md_log.h>    // md_log_type_t

#include <stdarg.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

struct md_allocator_i;

/*

WRITERS

Where output goes. A procedure which produces text (or bytes) writes it once, to a writer interface (md_writer_i), and
its caller decides where it ends up by which writer it hands over:

    md_str_writer_t     a string in memory, grown with an allocator
    md_buf_writer_t     the caller's fixed array, with snprintf semantics
    md_file_writer_t    a file, through a staging buffer
    md_log_writer_t     the log, a message per line
    md_writer_stdout()  stdout and stderr, written through
    md_writer_stderr()

    md_str_writer_t sw = md_str_writer_create(alloc);
    md_str_writer_fmt(&sw, "%d atoms\n", n);        // the string writer's own API
    md_script_ir_write_diagnostics(&sw.w, ir);      // handed over as the interface
    str_t s = md_str_writer_to_str(&sw);            // always current: there is nothing to flush
    ...
    md_str_writer_free(&sw);

    md_script_ir_write_diagnostics(md_writer_stdout(), ir);

THE INTERFACE

Every writer has an md_writer_i as its first member, w, and &x.w is what is handed over. The interface holds a buffer
[buf, end[ and a cursor, cur, in it: writing is a bounds check and a copy, inline, into the memory of the writer itself.
Only a write which does not fit consults the kind of writer, through its vtable: a string grows, a file passes on what
it has buffered, a fixed array drops what does not fit.

Through the interface a writer can only be written to: appended to, formatted into, flushed, asked how much has been
written. What only some writers can do (read back what is written, take it back, close a file) is done through the
concrete type, by whoever made it.

A procedure which is given an md_writer_i* only writes to it. It does not keep it once it returns, and if it opens temp
scopes while it writes, it opens them with md_temp_begin_avoid(md_writer_allocator(w)), as it would for an allocator it
is given: a string writer can grow in a temp arena, and a scope opened in the same arena would rewind over what it grows
into.

BUFFERS

Where the buffer comes from is up to each writer, and makes no difference to how fast writing is:
  - A string writer's buffer is the string itself.
  - A fixed array writer's buffer is the caller's array.
  - A file or log writer stages through a buffer: the caller's (a stack array, say), or one it allocates and frees.
  - The stdout and stderr writers have none. They write through to the C runtime's stdout and stderr, which buffer
    already, so what they write is ordered with printf, and is not lost if the program stops right after.

A writer is used by one thread at a time. None holds a pointer into itself, so they can be returned and copied by
value, but a copy is a second cursor into the same buffer: use one or the other.

*/

// The staging buffer a file or log writer allocates, when it is not given one
#define MD_WRITER_DEFAULT_CAPACITY (64 * 1024)

#if defined(__GNUC__) || defined(__clang__)
#define MD_WRITER_PRINTF(fmt_idx, arg_idx) __attribute__((format(printf, fmt_idx, arg_idx)))
#else
#define MD_WRITER_PRINTF(fmt_idx, arg_idx)
#endif

typedef struct md_writer_i md_writer_i;

typedef enum md_writer_flags_t {
    MD_WRITER_FLAG_ERROR     = 1u << 0,     // The sink failed. It sticks: what is written afterwards is discarded
    MD_WRITER_FLAG_TRUNCATED = 1u << 1,     // A fixed array writer dropped what did not fit
    MD_WRITER_FLAG_TTY       = 1u << 2,     // It writes to a terminal, where colour can be used
} md_writer_flags_t;

// What makes a kind of writer: the two slow paths. Implemented for each writer, and for a writer of one's own.
typedef struct md_writer_vtable_t {
    // Consumes what is buffered, [buf, cur[, followed by data[0, len[: passes it on, keeps it (growing the buffer) or
    // drops it, adds what it passed on or dropped to drained, and leaves room in [cur, end[ again.
    // With len 0 it flushes. Returns false if the sink failed (what it was given is lost).
    bool (*drain)(md_writer_i* w, const char* data, size_t len);

    // Makes room for n contiguous bytes at cur if it can, by draining or by growing. Afterwards [cur, end[ may still
    // be shorter than n, for a sink of fixed size. Returns false if the sink failed. NULL: it never makes room.
    bool (*rebase)(md_writer_i* w, size_t n);

    const char* name;   // "string", "file", ...: for the debugger
} md_writer_vtable_t;

struct md_writer_i {
    char* buf;          // The buffer [buf, end[. Never NULL: a writer without a buffer points at an empty one
    char* cur;          // [buf, cur[ is written and not yet passed on
    char* end;

    const md_writer_vtable_t* vt;
    struct md_allocator_i* alloc;   // What the writer allocated buf with (and grows it with), NULL if it was given buf
    size_t   drained;               // What was consumed before buf: drained + (cur - buf) bytes are written in all
    uint32_t flags;                 // md_writer_flags_t
};

// ### CONCRETE WRITERS ###

typedef struct md_str_writer_t {
    md_writer_i w;
#ifdef __cplusplus
    // Appending, and reading what is written as a str_t
    md_str_writer_t& operator+=(str_t str);
    md_str_writer_t& operator+=(const char* cstr);
    md_str_writer_t& operator+=(char c);
    operator str_t() const;
    operator md_writer_i&() { return w; }
#endif
} md_str_writer_t;

typedef struct md_buf_writer_t {
    md_writer_i w;
} md_buf_writer_t;

typedef struct md_file_writer_t {
    md_writer_i w;
    md_file_t file;
} md_file_writer_t;

typedef struct md_log_writer_t {
    md_writer_i w;
    md_log_type_t type;
} md_log_writer_t;

#ifdef __cplusplus
extern "C" {
#endif

// ### THE INTERFACE ###

// For a writer of one's own: a struct with an md_writer_i as its first member, whose vtable casts w back to it.
// buf and cap can be NULL and 0 for a writer without a buffer, whose writes all go to drain.
void md_writer_init(md_writer_i* w, const md_writer_vtable_t* vt, char* buf, size_t cap);

// The parts of md_write and md_writer_reserve which do not fit in the buffer. Not to be called directly.
void  md_writer_write_slow(md_writer_i* w, const void* data, size_t len);
char* md_writer_reserve_slow(md_writer_i* w, size_t n);

static inline void md_write(md_writer_i* w, const void* data, size_t len) {
    ASSERT(w);
    if (len <= (size_t)(w->end - w->cur)) {
        if (len) MEMCPY(w->cur, data, len);
        w->cur += len;
    } else {
        md_writer_write_slow(w, data, len);
    }
}

static inline void md_write_char(md_writer_i* w, char c) {
    ASSERT(w);
    if (w->cur < w->end) {
        *w->cur++ = c;
    } else {
        md_writer_write_slow(w, &c, 1);
    }
}

static inline void md_write_str(md_writer_i* w, str_t str) {
    md_write(w, str.ptr, str.len);
}

static inline void md_write_cstr(md_writer_i* w, const char* cstr) {
    if (cstr) md_write(w, cstr, strlen(cstr));
}

// c, count times
void md_write_fill(md_writer_i* w, char c, size_t count);

// Formatted as printf does. It formats in place when it fits, so a format is not measured first.
void md_writef (md_writer_i* w, const char* format, ...) MD_WRITER_PRINTF(2, 3);
void md_vwritef(md_writer_i* w, const char* format, va_list args);

// Integers, formatted without printf: for writing many of them
void md_write_i64(md_writer_i* w, int64_t x);
void md_write_u64(md_writer_i* w, uint64_t x);

// In fixed notation with decimals digits after the point, exactly as "%.*f" writes it
void md_write_f64(md_writer_i* w, double x, int decimals);

// Room for n contiguous bytes at the cursor, to format into in place: write at most n bytes there, then commit how
// many. NULL if the writer cannot make that much room (n is more than its buffer holds, or it has none): write the
// same through md_write then.
static inline char* md_writer_reserve(md_writer_i* w, size_t n) {
    ASSERT(w);
    return n <= (size_t)(w->end - w->cur) ? w->cur : md_writer_reserve_slow(w, n);
}

static inline void md_writer_commit(md_writer_i* w, size_t n) {
    ASSERT(w);
    ASSERT(n <= (size_t)(w->end - w->cur));
    w->cur += n;
}

// Passes on what is buffered: to the file, the log, or the C runtime (which is flushed as well). A writer in memory
// has nothing to pass on.
void md_writer_flush(md_writer_i* w);

// How many bytes have been written in all. For a fixed array writer this counts what was dropped as well.
static inline size_t md_writer_tell(const md_writer_i* w) {
    ASSERT(w);
    return w->drained + (size_t)(w->cur - w->buf);
}

// False once the sink has failed
static inline bool md_writer_ok(const md_writer_i* w) {
    ASSERT(w);
    return !(w->flags & MD_WRITER_FLAG_ERROR);
}

// The allocator the writer allocated its buffer with, NULL if it was given one: for md_temp_begin_avoid (see above)
static inline struct md_allocator_i* md_writer_allocator(const md_writer_i* w) {
    ASSERT(w);
    return w->alloc;
}

// Write through to the C runtime's stdout and stderr: no buffer, nothing to close, ordered with printf.
// There is one of each per thread, flagged MD_WRITER_FLAG_TTY when it writes to a terminal.
md_writer_i* md_writer_stdout(void);
md_writer_i* md_writer_stderr(void);

// ### STRING WRITER ###
// A string in memory, grown with an allocator, which can be a temp arena (see above).
// What is written can be read at any time, and is terminated when it is read.

void md_str_writer_init(md_str_writer_t* sw, struct md_allocator_i* alloc);
void md_str_writer_free(md_str_writer_t* sw);

static inline md_str_writer_t md_str_writer_create(struct md_allocator_i* alloc) {
    md_str_writer_t sw;
    md_str_writer_init(&sw, alloc);
    return sw;
}

static inline void md_str_writer_push_char (md_str_writer_t* sw, char c)                     { md_write_char(&sw->w, c); }
static inline void md_str_writer_push_str  (md_str_writer_t* sw, str_t str)                  { md_write_str(&sw->w, str); }
static inline void md_str_writer_push_cstr (md_str_writer_t* sw, const char* cstr)           { md_write_cstr(&sw->w, cstr); }
static inline void md_str_writer_push_cstrl(md_str_writer_t* sw, const char* cstr, size_t n) { md_write(&sw->w, cstr, n); }
void md_str_writer_fmt (md_str_writer_t* sw, const char* format, ...) MD_WRITER_PRINTF(2, 3);
void md_str_writer_vfmt(md_str_writer_t* sw, const char* format, va_list args);

static inline size_t md_str_writer_len(const md_str_writer_t* sw) {
    return (size_t)(sw->w.cur - sw->w.buf);
}

static inline bool md_str_writer_empty(const md_str_writer_t* sw) {
    return sw->w.cur == sw->w.buf;
}

// What is written, terminated (in a byte kept for it after what is written), so that its ptr is a C string as well.
// Valid until the next write, which can move the buffer as it grows.
static inline str_t md_str_writer_to_str(const md_str_writer_t* sw) {
    if (sw->w.end == sw->w.buf) {
        // Nothing allocated yet (a buffer, once allocated, always has room beyond the terminator)
        str_t empty = { "", 0 };
        return empty;
    }
    *sw->w.cur = '\0';
    str_t str = { sw->w.buf, (size_t)(sw->w.cur - sw->w.buf) };
    return str;
}

static inline const char* md_str_writer_to_cstr(const md_str_writer_t* sw) {
    return md_str_writer_to_str(sw).ptr;
}

// Takes back the last n bytes (all of them, if there are fewer)
void md_str_writer_pop(md_str_writer_t* sw, size_t n);

// Empties it. The buffer is kept.
void md_str_writer_reset(md_str_writer_t* sw);

// Hands over what is written as a string of its own: allocated with the writer's allocator, terminated, and freed with
// str_free. The writer is left empty and can be written to again, or freed (which then has nothing to free).
str_t md_str_writer_release(md_str_writer_t* sw);

// ### FIXED ARRAY WRITER ###
// The caller's array buf of cap bytes (at least 1), of which one is kept for the terminator. What does not fit is
// dropped, and flags MD_WRITER_FLAG_TRUNCATED, but is still counted: md_writer_tell gives the size it would have
// needed, as snprintf does. Nothing to free.

void md_buf_writer_init(md_buf_writer_t* bw, char* buf, size_t cap);

static inline md_buf_writer_t md_buf_writer_create(char* buf, size_t cap) {
    md_buf_writer_t bw;
    md_buf_writer_init(&bw, buf, cap);
    return bw;
}

// What is kept, terminated in the array
str_t md_buf_writer_to_str(const md_buf_writer_t* bw);

static inline bool md_buf_writer_truncated(const md_buf_writer_t* bw) {
    return (bw->w.flags & MD_WRITER_FLAG_TRUNCATED) != 0;
}

void md_buf_writer_reset(md_buf_writer_t* bw);

// ### FILE WRITER ###
// A file, through a staging buffer: the caller's buf of cap bytes, or with buf NULL, one of cap bytes
// (MD_WRITER_DEFAULT_CAPACITY if 0) which the writer allocates. The file stays the caller's: closing the writer flushes
// and frees the buffer, and does not close the file.

void md_file_writer_init(md_file_writer_t* fw, md_file_t file, char* buf, size_t cap);

static inline md_file_writer_t md_file_writer_create(md_file_t file, char* buf, size_t cap) {
    md_file_writer_t fw;
    md_file_writer_init(&fw, file, buf, cap);
    return fw;
}

// Flushes and frees. False if anything written failed to reach the file.
bool md_file_writer_close(md_file_writer_t* fw);

// ### LOG WRITER ###
// The log: each complete line is logged as a message of type, without its newline. buf and cap as for a file writer;
// a line longer than the buffer is logged in pieces, and a flush logs what there is of the last line.

void md_log_writer_init(md_log_writer_t* lw, md_log_type_t type, char* buf, size_t cap);

static inline md_log_writer_t md_log_writer_create(md_log_type_t type, char* buf, size_t cap) {
    md_log_writer_t lw;
    md_log_writer_init(&lw, type, buf, cap);
    return lw;
}

// Logs what remains and frees
void md_log_writer_close(md_log_writer_t* lw);

#ifdef __cplusplus
}

// Appending with operators: on the interface, and through it on any writer which converts to it.
// The integer overloads are exact matches for each type, so that w << 42 is not ambiguous; a char is a character.
inline md_writer_i& operator<<(md_writer_i& w, str_t str)            { md_write_str(&w, str);   return w; }
inline md_writer_i& operator<<(md_writer_i& w, const char* cstr)     { md_write_cstr(&w, cstr); return w; }
inline md_writer_i& operator<<(md_writer_i& w, char c)               { md_write_char(&w, c);    return w; }
inline md_writer_i& operator<<(md_writer_i& w, int x)                { md_write_i64(&w, x);     return w; }
inline md_writer_i& operator<<(md_writer_i& w, long x)               { md_write_i64(&w, x);     return w; }
inline md_writer_i& operator<<(md_writer_i& w, long long x)          { md_write_i64(&w, x);     return w; }
inline md_writer_i& operator<<(md_writer_i& w, unsigned x)           { md_write_u64(&w, x);     return w; }
inline md_writer_i& operator<<(md_writer_i& w, unsigned long x)      { md_write_u64(&w, x);     return w; }
inline md_writer_i& operator<<(md_writer_i& w, unsigned long long x) { md_write_u64(&w, x);     return w; }

inline md_str_writer_t& md_str_writer_t::operator+=(str_t str)        { md_write_str(&w, str);   return *this; }
inline md_str_writer_t& md_str_writer_t::operator+=(const char* cstr) { md_write_cstr(&w, cstr); return *this; }
inline md_str_writer_t& md_str_writer_t::operator+=(char c)           { md_write_char(&w, c);    return *this; }
inline md_str_writer_t::operator str_t() const                        { return md_str_writer_to_str(this); }
#endif
