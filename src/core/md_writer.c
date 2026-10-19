#include <core/md_writer.h>

#include <core/md_allocator.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_platform.h>

#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#if MD_PLATFORM_WINDOWS
#include <io.h>
#define md_isatty(f) (_isatty(_fileno(f)) != 0)
#else
#include <unistd.h>
#define md_isatty(f) (isatty(fileno(f)) != 0)
#endif

// The writers are handed over as &x.w, and their vtables cast w back to the writer
STATIC_ASSERT(offsetof(md_str_writer_t,  w) == 0, "md_writer_i must be the first member");
STATIC_ASSERT(offsetof(md_buf_writer_t,  w) == 0, "md_writer_i must be the first member");
STATIC_ASSERT(offsetof(md_file_writer_t, w) == 0, "md_writer_i must be the first member");
STATIC_ASSERT(offsetof(md_log_writer_t,  w) == 0, "md_writer_i must be the first member");

// What a writer without a buffer points at. It has no room, so it is never written to.
static char empty_buf[1];

static inline size_t room(const md_writer_i* w) {
    return (size_t)(w->end - w->cur);
}

static inline bool has_buffer(const md_writer_i* w) {
    return w->buf && w->buf != empty_buf;
}

// ### THE INTERFACE ###

void md_writer_init(md_writer_i* w, const md_writer_vtable_t* vt, char* buf, size_t cap) {
    ASSERT(w);
    ASSERT(vt && vt->drain);
    MEMSET(w, 0, sizeof(md_writer_i));
    w->vt = vt;
    if (buf && cap) {
        w->buf = w->cur = buf;
        w->end = buf + cap;
    } else {
        w->buf = w->cur = w->end = empty_buf;
    }
}

// Once the sink has failed, what is written is discarded, but still counted
static void discard(md_writer_i* w, size_t len) {
    w->drained += (size_t)(w->cur - w->buf) + len;
    w->cur = w->buf;
}

void md_writer_write_slow(md_writer_i* w, const void* data, size_t len) {
    ASSERT(w && w->vt && w->vt->drain);
    if (w->flags & MD_WRITER_FLAG_ERROR) {
        discard(w, len);
        return;
    }
    if (!w->vt->drain(w, (const char*)data, len)) {
        w->flags |= MD_WRITER_FLAG_ERROR;
    }
}

char* md_writer_reserve_slow(md_writer_i* w, size_t n) {
    ASSERT(w && w->vt);
    if (w->flags & MD_WRITER_FLAG_ERROR) {
        discard(w, 0);
    } else if (w->vt->rebase && !w->vt->rebase(w, n)) {
        w->flags |= MD_WRITER_FLAG_ERROR;
        discard(w, 0);
    }
    return n <= room(w) ? w->cur : NULL;
}

void md_writer_flush(md_writer_i* w) {
    ASSERT(w && w->vt && w->vt->drain);
    if (w->flags & MD_WRITER_FLAG_ERROR) {
        discard(w, 0);
        return;
    }
    if (!w->vt->drain(w, NULL, 0)) {
        w->flags |= MD_WRITER_FLAG_ERROR;
    }
}

void md_write_fill(md_writer_i* w, char c, size_t count) {
    ASSERT(w);
    while (count > 0) {
        const size_t n = MIN(room(w), count);
        if (n > 0) {
            MEMSET(w->cur, c, n);
            w->cur += n;
            count -= n;
        } else {
            char block[64];
            const size_t m = MIN(sizeof(block), count);
            MEMSET(block, c, m);
            md_writer_write_slow(w, block, m);
            count -= m;
        }
    }
}

void md_vwritef(md_writer_i* w, const char* format, va_list args) {
    ASSERT(w);
    ASSERT(format);

    // Formatted straight into the buffer: most formats fit
    const size_t avail = room(w);
    va_list probe;
    va_copy(probe, args);
    const int res = md_vsnprintf(w->cur, avail, format, probe);
    va_end(probe);

    if (res < 0) {
        MD_LOG_ERROR("Formatting failed: '%s'", format);
        return;
    }

    const size_t len = (size_t)res;
    if (len < avail) {
        w->cur += len;
        return;
    }

    // It did not fit, terminator included: make room for it and format it again
    char* dst = md_writer_reserve(w, len + 1);
    if (dst) {
        md_vsnprintf(dst, len + 1, format, args);
        w->cur += len;
        return;
    }

    // The writer cannot hold it in one piece (a fixed buffer, or none): format it elsewhere and write it through
    char stack_buf[1024];
    if (len < sizeof(stack_buf)) {
        md_vsnprintf(stack_buf, sizeof(stack_buf), format, args);
        md_write(w, stack_buf, len);
    } else {
        md_allocator_i* heap = md_get_heap_allocator();
        char* tmp = (char*)md_alloc(heap, len + 1);
        md_vsnprintf(tmp, len + 1, format, args);
        md_write(w, tmp, len);
        md_free(heap, tmp, len + 1);
    }
}

void md_writef(md_writer_i* w, const char* format, ...) {
    va_list args;
    va_start(args, format);
    md_vwritef(w, format, args);
    va_end(args);
}

void md_write_u64(md_writer_i* w, uint64_t x) {
    char tmp[20];   // UINT64_MAX has 20 digits
    char* p = tmp + sizeof(tmp);
    do {
        *--p = (char)('0' + x % 10);
        x /= 10;
    } while (x);
    md_write(w, p, (size_t)(tmp + sizeof(tmp) - p));
}

void md_write_i64(md_writer_i* w, int64_t x) {
    char tmp[21];   // INT64_MIN has 19 digits and a sign
    char* p = tmp + sizeof(tmp);
    uint64_t u = x < 0 ? (uint64_t)0 - (uint64_t)x : (uint64_t)x;
    do {
        *--p = (char)('0' + u % 10);
        u /= 10;
    } while (u);
    if (x < 0) *--p = '-';
    md_write(w, p, (size_t)(tmp + sizeof(tmp) - p));
}

void md_write_f64(md_writer_i* w, double x, int decimals) {
    // Through printf for now: it rounds correctly, which a faster version has to as well before it replaces this
    md_writef(w, "%.*f", MAX(decimals, 0), x);
}

// ### STDOUT AND STDERR ###
// No buffer: every write is data for drain, which writes it through to the C runtime

static bool std_drain(md_writer_i* w, FILE* stream, const char* data, size_t len) {
    if (len == 0) {
        return fflush(stream) == 0;
    }
    w->drained += len;
    return fwrite(data, 1, len, stream) == len;
}

static bool stdout_drain(md_writer_i* w, const char* data, size_t len) { return std_drain(w, stdout, data, len); }
static bool stderr_drain(md_writer_i* w, const char* data, size_t len) { return std_drain(w, stderr, data, len); }

static const md_writer_vtable_t stdout_vt = { stdout_drain, NULL, "stdout" };
static const md_writer_vtable_t stderr_vt = { stderr_drain, NULL, "stderr" };

static THREAD_LOCAL md_writer_i tl_stdout_writer;
static THREAD_LOCAL md_writer_i tl_stderr_writer;

md_writer_i* md_writer_stdout(void) {
    md_writer_i* w = &tl_stdout_writer;
    if (!w->vt) {
        md_writer_init(w, &stdout_vt, NULL, 0);
        if (md_isatty(stdout)) w->flags |= MD_WRITER_FLAG_TTY;
    }
    return w;
}

md_writer_i* md_writer_stderr(void) {
    md_writer_i* w = &tl_stderr_writer;
    if (!w->vt) {
        md_writer_init(w, &stderr_vt, NULL, 0);
        if (md_isatty(stderr)) w->flags |= MD_WRITER_FLAG_TTY;
    }
    return w;
}

// ### STRING WRITER ###
// The allocation is [buf, end] (end included): the byte at end is kept for the terminator.

static size_t str_alloc_size(const md_writer_i* w) {
    return has_buffer(w) ? (size_t)(w->end - w->buf) + 1 : 0;
}

// Room for n more bytes at cur
static bool str_grow(md_writer_i* w, size_t n) {
    ASSERT(w->alloc);
    const size_t len = (size_t)(w->cur - w->buf);
    const size_t old_size = str_alloc_size(w);
    size_t new_size = MAX(old_size * 2, 64);
    if (new_size < len + n + 1) {
        new_size = len + n + 1;
    }
    char* buf = old_size ? (char*)md_realloc(w->alloc, w->buf, old_size, new_size) : (char*)md_alloc(w->alloc, new_size);
    if (!buf) {
        return false;
    }
    w->buf = buf;
    w->cur = buf + len;
    w->end = buf + new_size - 1;
    return true;
}

static bool str_drain(md_writer_i* w, const char* data, size_t len) {
    // Nothing is passed on: it is all kept
    if (len == 0) {
        return true;
    }
    if (room(w) < len && !str_grow(w, len)) {
        return false;
    }
    MEMCPY(w->cur, data, len);
    w->cur += len;
    return true;
}

static bool str_rebase(md_writer_i* w, size_t n) {
    return room(w) >= n || str_grow(w, n);
}

static const md_writer_vtable_t str_vt = { str_drain, str_rebase, "string" };

void md_str_writer_init(md_str_writer_t* sw, md_allocator_i* alloc) {
    ASSERT(sw);
    ASSERT(alloc);
    md_writer_init(&sw->w, &str_vt, NULL, 0);
    sw->w.alloc = alloc;
}

void md_str_writer_free(md_str_writer_t* sw) {
    ASSERT(sw);
    if (has_buffer(&sw->w)) {
        ASSERT(sw->w.alloc);
        md_free(sw->w.alloc, sw->w.buf, str_alloc_size(&sw->w));
    }
    MEMSET(sw, 0, sizeof(md_str_writer_t));
}

void md_str_writer_vfmt(md_str_writer_t* sw, const char* format, va_list args) {
    ASSERT(sw);
    md_vwritef(&sw->w, format, args);
}

void md_str_writer_fmt(md_str_writer_t* sw, const char* format, ...) {
    va_list args;
    va_start(args, format);
    md_vwritef(&sw->w, format, args);
    va_end(args);
}

void md_str_writer_pop(md_str_writer_t* sw, size_t n) {
    ASSERT(sw);
    sw->w.cur -= MIN(n, md_str_writer_len(sw));
}

void md_str_writer_reset(md_str_writer_t* sw) {
    ASSERT(sw);
    sw->w.cur = sw->w.buf;
    sw->w.drained = 0;
    sw->w.flags = 0;
}

str_t md_str_writer_release(md_str_writer_t* sw) {
    ASSERT(sw);
    str_t str = {0};
    md_allocator_i* alloc = sw->w.alloc;
    if (has_buffer(&sw->w)) {
        const size_t len  = md_str_writer_len(sw);
        const size_t size = str_alloc_size(&sw->w);
        if (len > 0) {
            // Down to what str_free frees: the string and its terminator
            char* ptr = (char*)md_realloc(alloc, sw->w.buf, size, len + 1);
            ptr[len] = '\0';
            str.ptr = ptr;
            str.len = len;
        } else {
            md_free(alloc, sw->w.buf, size);
        }
    }
    md_str_writer_init(sw, alloc);
    return str;
}

// ### FIXED ARRAY WRITER ###
// The byte at end is kept for the terminator, as for a string writer

static bool buf_drain(md_writer_i* w, const char* data, size_t len) {
    // Keeps what fits. The rest is dropped, but counted
    const size_t n = MIN(room(w), len);
    if (n > 0) {
        MEMCPY(w->cur, data, n);
        w->cur += n;
    }
    if (n < len) {
        w->drained += len - n;
        w->flags |= MD_WRITER_FLAG_TRUNCATED;
    }
    return true;
}

// No rebase: there is no more room to be made
static const md_writer_vtable_t buf_vt = { buf_drain, NULL, "array" };

void md_buf_writer_init(md_buf_writer_t* bw, char* buf, size_t cap) {
    ASSERT(bw);
    ASSERT(buf);
    ASSERT(cap >= 1);
    md_writer_init(&bw->w, &buf_vt, NULL, 0);
    bw->w.buf = bw->w.cur = buf;
    bw->w.end = buf + cap - 1;
    *buf = '\0';
}

str_t md_buf_writer_to_str(const md_buf_writer_t* bw) {
    ASSERT(bw);
    *bw->w.cur = '\0';
    str_t str = { bw->w.buf, (size_t)(bw->w.cur - bw->w.buf) };
    return str;
}

void md_buf_writer_reset(md_buf_writer_t* bw) {
    ASSERT(bw);
    bw->w.cur = bw->w.buf;
    bw->w.drained = 0;
    bw->w.flags = 0;
    *bw->w.buf = '\0';
}

// ### FILE WRITER ###

static bool file_drain(md_writer_i* w, const char* data, size_t len) {
    md_file_writer_t* fw = (md_file_writer_t*)w;
    bool ok = true;

    const size_t buffered = (size_t)(w->cur - w->buf);
    if (buffered > 0) {
        ok = md_file_write(fw->file, w->buf, buffered) == buffered;
        w->drained += buffered;
        w->cur = w->buf;
    }

    if (len > 0) {
        if (len < (size_t)(w->end - w->buf)) {
            // Small enough to be buffered until the next drain
            MEMCPY(w->cur, data, len);
            w->cur += len;
        } else {
            ok = (md_file_write(fw->file, data, len) == len) && ok;
            w->drained += len;
        }
    }
    return ok;
}

static bool file_rebase(md_writer_i* w, size_t n) {
    (void)n;
    return file_drain(w, NULL, 0);
}

static const md_writer_vtable_t file_vt = { file_drain, file_rebase, "file" };

void md_file_writer_init(md_file_writer_t* fw, md_file_t file, char* buf, size_t cap) {
    ASSERT(fw);
    ASSERT(md_file_valid(file));
    md_allocator_i* alloc = NULL;
    if (!buf) {
        cap = cap ? cap : MD_WRITER_DEFAULT_CAPACITY;
        alloc = md_get_heap_allocator();
        buf = (char*)md_alloc(alloc, cap);
    }
    md_writer_init(&fw->w, &file_vt, buf, cap);
    fw->w.alloc = alloc;
    fw->file = file;
}

bool md_file_writer_close(md_file_writer_t* fw) {
    ASSERT(fw);
    md_writer_flush(&fw->w);
    const bool ok = md_writer_ok(&fw->w);
    if (fw->w.alloc && has_buffer(&fw->w)) {
        md_free(fw->w.alloc, fw->w.buf, (size_t)(fw->w.end - fw->w.buf));
    }
    MEMSET(fw, 0, sizeof(md_file_writer_t));
    return ok;
}

// ### LOG WRITER ###
// The byte at end is kept for terminating the last line, when it is logged unfinished

// Logs every complete line in [buf, cur[, and with all, what there is after the last one too. Keeps the rest.
static void log_lines(md_log_writer_t* lw, bool all) {
    md_writer_i* w = &lw->w;
    char* line = w->buf;
    for (char* p = w->buf; p < w->cur; ++p) {
        if (*p == '\n') {
            char* e = (p > line && p[-1] == '\r') ? p - 1 : p;
            *e = '\0';
            md_log(lw->type, line);
            line = p + 1;
        }
    }
    if (all && line < w->cur) {
        *w->cur = '\0';
        md_log(lw->type, line);
        line = w->cur;
    }
    const size_t consumed = (size_t)(line - w->buf);
    if (consumed > 0) {
        const size_t rest = (size_t)(w->cur - line);
        MEMMOVE(w->buf, line, rest);
        w->cur = w->buf + rest;
        w->drained += consumed;
    }
}

static bool log_drain(md_writer_i* w, const char* data, size_t len) {
    md_log_writer_t* lw = (md_log_writer_t*)w;
    // A flush logs what there is of the last line too
    log_lines(lw, len == 0);
    while (len > 0) {
        if (w->cur == w->end) {
            // A line longer than the buffer: logged in pieces
            log_lines(lw, true);
        }
        const size_t n = MIN(room(w), len);
        MEMCPY(w->cur, data, n);
        w->cur += n;
        data += n;
        len -= n;
        log_lines(lw, false);
    }
    return true;
}

static bool log_rebase(md_writer_i* w, size_t n) {
    (void)n;
    log_lines((md_log_writer_t*)w, false);
    return true;
}

static const md_writer_vtable_t log_vt = { log_drain, log_rebase, "log" };

void md_log_writer_init(md_log_writer_t* lw, md_log_type_t type, char* buf, size_t cap) {
    ASSERT(lw);
    md_allocator_i* alloc = NULL;
    if (!buf) {
        cap = cap ? cap : MD_WRITER_DEFAULT_CAPACITY;
        alloc = md_get_heap_allocator();
        buf = (char*)md_alloc(alloc, cap);
    }
    ASSERT(cap >= 2);   // a byte of a line, and its terminator
    md_writer_init(&lw->w, &log_vt, buf, cap);
    lw->w.end -= 1;
    lw->w.alloc = alloc;
    lw->type = type;
}

void md_log_writer_close(md_log_writer_t* lw) {
    ASSERT(lw);
    md_writer_flush(&lw->w);
    if (lw->w.alloc && has_buffer(&lw->w)) {
        md_free(lw->w.alloc, lw->w.buf, (size_t)(lw->w.end - lw->w.buf) + 1);
    }
    MEMSET(lw, 0, sizeof(md_log_writer_t));
}
