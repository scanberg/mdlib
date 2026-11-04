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

// Keeps what is not taken at the start of the buffer, and reads as much after it as fits. Returns the bytes read: 0 at
// the end of the source, for a string, or when the buffer is full of what is not taken.
static size_t refill(md_reader_i* r) {
    if (!r->vt || (r->flags & MD_READER_FLAG_END)) return 0;
    const size_t rest = (size_t)(r->end - r->cur);
    if (rest == r->cap) return 0;
    r->base += (int64_t)(r->cur - r->buf);
    if (rest && r->cur != r->mem) memmove(r->mem, r->cur, rest);
    r->buf = r->cur = r->mem;
    r->end = r->mem + rest;

    const size_t want = r->cap - rest;
    const size_t n = r->vt->read(r, r->mem + rest, want);
    r->end += n;
    if (n < want) r->flags |= MD_READER_FLAG_END;
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
        const char* nl = (const char*)memchr(from, '\n', (size_t)(r->end - from));
        if (nl) {
            *line = md_reader_strip_cr((str_t){ r->cur, (size_t)(nl - r->cur) });
            *next = nl + 1;
            return true;
        }
        searched = (size_t)(r->end - r->cur);
    }
}

str_t md_reader_block(md_reader_i* r) {
    ASSERT(r);
    refill(r);
    const char* text = r->cur;
    size_t len = (size_t)(r->end - r->cur);
    if (!(r->flags & MD_READER_FLAG_END)) {
        // Whole lines: up to the last newline (the rest is read again with the next block)
        const char* c = text + len;
        while (c > text && c[-1] != '\n') --c;
        if (c > text) len = (size_t)(c - text);     // else a line longer than the buffer: cut at it
    }
    return (str_t){ text, len };
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

static size_t file_read(md_reader_i* r, char* dst, size_t cap) {
    md_file_reader_t* fr = (md_file_reader_t*)r;
    size_t total = 0;
    // A read may come back short of what was asked for before the end (large requests on some platforms)
    while (total < cap) {
        const size_t n = md_file_read(fr->file, dst + total, cap - total);
        if (n == 0) break;
        total += n;
    }
    return total;
}

static bool file_seek(md_reader_i* r, int64_t offset) {
    md_file_reader_t* fr = (md_file_reader_t*)r;
    return md_file_seek(fr->file, offset, MD_FILE_BEG);
}

static const md_reader_vtable_t file_vtable = {
    .read = file_read,
    .seek = file_seek,
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
    // From where the file is positioned
    fr->r.base = (int64_t)md_file_tell(file);
    return true;
}

bool md_file_reader_open(md_file_reader_t* fr, str_t path, size_t cap, struct md_allocator_i* alloc) {
    ASSERT(fr);
    md_file_t file = {0};
    if (!md_file_open(&file, path, MD_FILE_READ)) {
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
