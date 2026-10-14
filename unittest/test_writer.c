#include "utest.h"

#include <core/md_writer.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_str.h>

#include <stdio.h>
#include <string.h>

// Writes through the interface only, as a procedure which is handed a writer does
static void produce(md_writer_i* w) {
    md_write_cstr(w, "atoms: ");
    md_write_i64(w, 42);
    md_write_char(w, '\n');
    md_writef(w, "%s=%.2f\n", "x", 1.5);
}

UTEST(writer, str_writer) {
    md_str_writer_t sw = md_str_writer_create(md_get_heap_allocator());
    EXPECT_TRUE(md_str_writer_empty(&sw));
    EXPECT_STREQ("", md_str_writer_to_cstr(&sw));
    EXPECT_EQ((size_t)0, md_str_writer_to_str(&sw).len);
    EXPECT_STREQ("", md_str_writer_to_str(&sw).ptr);     // to_str is terminated too, empty or not

    md_str_writer_push_cstr(&sw, "Hej");
    md_str_writer_fmt(&sw, " %d Cool", 2);
    md_str_writer_push_str(&sw, STR_LIT(" 4 School."));
    EXPECT_STREQ("Hej 2 Cool 4 School.", md_str_writer_to_cstr(&sw));
    EXPECT_STREQ("Hej 2 Cool 4 School.", md_str_writer_to_str(&sw).ptr);
    EXPECT_EQ((size_t)20, md_str_writer_len(&sw));

    // Written through the interface, read through the string writer: one cursor, nothing to flush
    produce(&sw.w);
    EXPECT_STREQ("Hej 2 Cool 4 School.atoms: 42\nx=1.50\n", md_str_writer_to_cstr(&sw));
    EXPECT_EQ(md_str_writer_len(&sw), md_writer_tell(&sw.w));
    EXPECT_TRUE(md_writer_ok(&sw.w));

    md_str_writer_pop(&sw, 1);
    str_t str = md_str_writer_to_str(&sw);
    EXPECT_TRUE(str_eq(str_substr(str, str.len - 6, 6), STR_LIT("x=1.50")));
    md_str_writer_pop(&sw, 1000);
    EXPECT_TRUE(md_str_writer_empty(&sw));

    md_str_writer_push_cstrl(&sw, "abcdef", 3);
    md_str_writer_push_char(&sw, '!');
    EXPECT_STREQ("abc!", md_str_writer_to_cstr(&sw));

    md_str_writer_reset(&sw);
    EXPECT_TRUE(md_str_writer_empty(&sw));
    EXPECT_STREQ("", md_str_writer_to_cstr(&sw));

    md_str_writer_free(&sw);
    md_str_writer_free(&sw);    // freeing twice, or a zeroed writer, is harmless

    md_str_writer_t zeroed = {0};
    EXPECT_STREQ("", md_str_writer_to_cstr(&zeroed));
    md_str_writer_free(&zeroed);
}

UTEST(writer, str_writer_growth) {
    md_allocator_i* heap = md_get_heap_allocator();
    md_str_writer_t sw = md_str_writer_create(heap);

    // Many small formats across many growths, most of which do not fit in what is left of the buffer
    const int count = 100000;
    for (int i = 0; i < count; ++i) {
        md_str_writer_fmt(&sw, "%d,", i);
    }
    str_t str = md_str_writer_to_str(&sw);
    int expected = 0;
    bool in_order = true;
    for (size_t i = 0; i < str.len;) {
        int value = 0;
        while (str.ptr[i] != ',') value = value * 10 + (str.ptr[i++] - '0');
        i += 1;
        in_order &= value == expected++;
    }
    EXPECT_TRUE(in_order);
    EXPECT_EQ(count, expected);

    // A format larger than anything buffered so far
    md_str_writer_reset(&sw);
    char big[5000];
    MEMSET(big, 'a', sizeof(big) - 1);
    big[sizeof(big) - 1] = '\0';
    md_str_writer_fmt(&sw, "<%s>", big);
    str = md_str_writer_to_str(&sw);
    EXPECT_EQ(sizeof(big) + 1, str.len);
    EXPECT_EQ('<', str.ptr[0]);
    EXPECT_EQ('a', str.ptr[str.len / 2]);
    EXPECT_EQ('>', str.ptr[str.len - 1]);

    md_str_writer_free(&sw);
}

UTEST(writer, str_writer_release) {
    md_allocator_i* heap  = md_get_heap_allocator();
    md_allocator_i* arena = md_arena_allocator_create(heap, KILOBYTES(64));

    md_str_writer_t sw = md_str_writer_create(arena);
    md_str_writer_push_cstr(&sw, "released");
    str_t released = md_str_writer_release(&sw);
    EXPECT_EQ((size_t)8, released.len);
    EXPECT_STREQ("released", released.ptr);    // terminated

    // The writer is left empty, and can be written to again
    EXPECT_TRUE(md_str_writer_empty(&sw));
    md_str_writer_push_cstr(&sw, "again");
    EXPECT_STREQ("again", md_str_writer_to_cstr(&sw));
    md_str_writer_free(&sw);
    EXPECT_STREQ("released", released.ptr);    // the released string is its own
    str_free(released, arena);

    md_str_writer_t hw = md_str_writer_create(heap);
    EXPECT_EQ((size_t)0, md_str_writer_release(&hw).len);
    md_str_writer_push_cstr(&hw, "x");
    str_t x = md_str_writer_release(&hw);
    EXPECT_TRUE(str_eq(x, STR_LIT("x")));
    str_free(x, heap);
    md_str_writer_free(&hw);

    md_arena_allocator_destroy(arena);
}

// Opens temp scopes while it writes, avoiding the writer's allocator as it should
static void produce_with_temp(md_writer_i* w, int count) {
    for (int i = 0; i < count; ++i) {
        md_temp_scope_t temp = md_temp_begin_avoid(md_writer_allocator(w));
        char* scratch = (char*)md_temp_alloc(temp, 256);
        snprintf(scratch, 256, "line %d", i);
        md_writef(w, "%s;", scratch);   // grows the writer, in the arena the scope avoids
        md_temp_end(temp);
    }
}

UTEST(writer, str_writer_temp_arena) {
    md_temp_scope_t temp = md_temp_begin();
    md_str_writer_t sw = md_str_writer_create(md_temp_allocator(temp));
    EXPECT_TRUE(md_writer_allocator(&sw.w) == md_temp_allocator(temp));

    md_str_writer_push_cstr(&sw, "start;");
    produce_with_temp(&sw.w, 200);

    md_str_writer_t expected = md_str_writer_create(md_get_heap_allocator());
    md_str_writer_push_cstr(&expected, "start;");
    for (int i = 0; i < 200; ++i) md_str_writer_fmt(&expected, "line %d;", i);
    EXPECT_TRUE(str_eq(md_str_writer_to_str(&expected), md_str_writer_to_str(&sw)));

    md_str_writer_free(&expected);
    md_temp_end(temp);
}

UTEST(writer, buf_writer) {
    char buf[8];
    md_buf_writer_t bw = md_buf_writer_create(buf, sizeof(buf));
    EXPECT_STREQ("", buf);

    md_write_cstr(&bw.w, "abc");
    EXPECT_FALSE(md_buf_writer_truncated(&bw));
    EXPECT_EQ((size_t)3, md_buf_writer_to_str(&bw).len);
    EXPECT_STREQ("abc", buf);

    // 7 of the 8 bytes hold text, one the terminator: of 9 digits, 4 fit
    md_writef(&bw.w, "%d", 123456789);
    EXPECT_TRUE(md_buf_writer_truncated(&bw));
    EXPECT_EQ((size_t)7, md_buf_writer_to_str(&bw).len);
    EXPECT_STREQ("abc1234", buf);
    EXPECT_EQ((size_t)12, md_writer_tell(&bw.w));     // the size it would have needed

    md_write_char(&bw.w, 'x');
    EXPECT_EQ((size_t)13, md_writer_tell(&bw.w));
    EXPECT_STREQ("abc1234", md_buf_writer_to_str(&bw).ptr);
    EXPECT_TRUE(md_writer_ok(&bw.w));                   // truncating is not failing

    md_buf_writer_reset(&bw);
    EXPECT_FALSE(md_buf_writer_truncated(&bw));
    EXPECT_EQ((size_t)0, md_writer_tell(&bw.w));

    // What fits exactly is not truncated
    md_writef(&bw.w, "%s", "1234567");
    EXPECT_FALSE(md_buf_writer_truncated(&bw));
    EXPECT_STREQ("1234567", md_buf_writer_to_str(&bw).ptr);
    EXPECT_TRUE(md_writer_reserve(&bw.w, 1) == NULL);   // full: no room can be made

    // A single byte holds only the terminator
    char one[1];
    md_buf_writer_t ow = md_buf_writer_create(one, sizeof(one));
    md_write_cstr(&ow.w, "x");
    EXPECT_TRUE(md_buf_writer_truncated(&ow));
    EXPECT_EQ((size_t)0, md_buf_writer_to_str(&ow).len);
    EXPECT_STREQ("", one);
}

UTEST(writer, numbers) {
    char buf[256];
    md_buf_writer_t bw = md_buf_writer_create(buf, sizeof(buf));
    md_write_i64(&bw.w, 0);         md_write_char(&bw.w, ' ');
    md_write_i64(&bw.w, -1);        md_write_char(&bw.w, ' ');
    md_write_i64(&bw.w, INT64_MIN); md_write_char(&bw.w, ' ');
    md_write_i64(&bw.w, INT64_MAX); md_write_char(&bw.w, ' ');
    md_write_u64(&bw.w, 0);         md_write_char(&bw.w, ' ');
    md_write_u64(&bw.w, UINT64_MAX);
    EXPECT_STREQ("0 -1 -9223372036854775808 9223372036854775807 0 18446744073709551615", md_buf_writer_to_str(&bw).ptr);

    // As printf writes them, whatever its rounding of halfway cases
    const double values[] = { 3.14159, -2.71828, 1.0e6, 0.0, -0.001 };
    const int decimals[]  = { 2, 3, 0, 4, 1 };
    for (size_t i = 0; i < ARRAY_SIZE(values); ++i) {
        char expected[64];
        snprintf(expected, sizeof(expected), "%.*f", decimals[i], values[i]);
        md_buf_writer_reset(&bw);
        md_write_f64(&bw.w, values[i], decimals[i]);
        EXPECT_STREQ(expected, md_buf_writer_to_str(&bw).ptr);
    }
}

UTEST(writer, fill_and_reserve) {
    md_str_writer_t sw = md_str_writer_create(md_get_heap_allocator());
    md_write_fill(&sw.w, '~', 1000);
    str_t str = md_str_writer_to_str(&sw);
    EXPECT_EQ((size_t)1000, str.len);
    EXPECT_EQ('~', str.ptr[0]);
    EXPECT_EQ('~', str.ptr[999]);

    char* p = md_writer_reserve(&sw.w, 4096);
    ASSERT_TRUE(p != NULL);
    MEMCPY(p, "0123456789", 10);
    md_writer_commit(&sw.w, 10);
    str = md_str_writer_to_str(&sw);
    EXPECT_EQ((size_t)1010, str.len);
    EXPECT_TRUE(str_eq(str_substr(str, 1000, 10), STR_LIT("0123456789")));
    md_str_writer_free(&sw);

    // A writer without a buffer cannot reserve
    EXPECT_TRUE(md_writer_reserve(md_writer_stdout(), 1) == NULL);
}

// Everything the file test writes, to whichever writer
static void produce_file_content(md_writer_i* w) {
    char block[100];
    MEMSET(block, '#', sizeof(block));
    for (int i = 0; i < 300; ++i) {
        if (i % 7 == 0) {
            md_writef(w, "%d:%s\n", i, "a line which is longer than the staging buffer of sixteen bytes");
        } else {
            md_writef(w, "%d\n", i);
        }
        if (i % 50 == 0) md_write(w, block, sizeof(block));   // written past the buffer
        if (i % 30 == 0) md_write_fill(w, '-', 40);
    }
}

UTEST(writer, file_writer) {
    md_allocator_i* heap = md_get_heap_allocator();
    const str_t path = STR_LIT("md_writer_test.tmp");

    md_str_writer_t expected = md_str_writer_create(heap);
    produce_file_content(&expected.w);

    // A staging buffer smaller than many of the writes, and one the writer allocates
    const size_t caps[] = { 16, 0 };
    char stage[16];
    for (size_t c = 0; c < ARRAY_SIZE(caps); ++c) {
        md_file_t file = {0};
        ASSERT_TRUE(md_file_open(&file, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE));
        md_file_writer_t fw = md_file_writer_create(file, caps[c] ? stage : NULL, caps[c]);
        produce_file_content(&fw.w);
        EXPECT_EQ(md_str_writer_len(&expected), md_writer_tell(&fw.w));
        EXPECT_TRUE(md_file_writer_close(&fw));
        md_file_close(&file);

        str_t written = load_textfile(path, heap);
        EXPECT_TRUE(str_eq(md_str_writer_to_str(&expected), written));
        str_free(written, heap);
    }

    md_str_writer_free(&expected);
    remove(path.ptr);
}

static struct {
    char lines[16][64];
    int count;
} log_capture;

static void capture_log(struct md_logger_o* inst, md_log_type_t type, const char* msg) {
    (void)inst;
    (void)type;
    if (log_capture.count < (int)ARRAY_SIZE(log_capture.lines)) {
        snprintf(log_capture.lines[log_capture.count++], sizeof(log_capture.lines[0]), "%s", msg);
    }
}

UTEST(writer, log_writer) {
    md_logger_i logger = { NULL, capture_log };
    MEMSET(&log_capture, 0, sizeof(log_capture));
    md_log_register(&logger);

    // 15 bytes of buffer: lines are logged as they complete, a longer one in pieces, the last one when closed
    char stage[16];
    md_log_writer_t lw = md_log_writer_create(MD_LOG_TYPE_DEBUG, stage, sizeof(stage));
    md_write_cstr(&lw.w, "first line\nsecond");
    md_write_cstr(&lw.w, " part\r\n");
    md_write_cstr(&lw.w, "a line longer than the staging buffer\n");
    md_write_cstr(&lw.w, "tail");
    md_log_writer_close(&lw);

    md_log_unregister(&logger);

    const char* expected[] = { "first line", "second part", "a line longer t", "han the staging", " buffer", "tail" };
    EXPECT_EQ((int)ARRAY_SIZE(expected), log_capture.count);
    for (int i = 0; i < (int)ARRAY_SIZE(expected) && i < log_capture.count; ++i) {
        EXPECT_STREQ(expected[i], log_capture.lines[i]);
    }
}

UTEST(writer, stdout_writer) {
    md_writer_i* out = md_writer_stdout();
    EXPECT_TRUE(out == md_writer_stdout());     // one per thread
    EXPECT_TRUE(out != md_writer_stderr());

    const char* line = "[ writer ] written through md_writer_stdout\n";
    const size_t before = md_writer_tell(out);
    md_writef(out, "%s", line);
    md_writer_flush(out);
    EXPECT_EQ(before + strlen(line), md_writer_tell(out));
    EXPECT_TRUE(md_writer_ok(out));
}

// A sink of one's own, which fails after a number of bytes
typedef struct failing_writer_t {
    md_writer_i w;
    size_t passed;
    size_t fail_after;
} failing_writer_t;

static bool failing_drain(md_writer_i* w, const char* data, size_t len) {
    (void)data;
    failing_writer_t* fw = (failing_writer_t*)w;
    const size_t n = (size_t)(w->cur - w->buf) + len;
    w->drained += n;
    w->cur = w->buf;
    if (fw->passed + n > fw->fail_after) return false;
    fw->passed += n;
    return true;
}

static const md_writer_vtable_t failing_vt = { failing_drain, NULL, "failing" };

UTEST(writer, custom_writer_and_errors) {
    char stage[8];
    failing_writer_t fw = { .fail_after = 20 };
    md_writer_init(&fw.w, &failing_vt, stage, sizeof(stage));

    md_write_cstr(&fw.w, "0123456789");
    EXPECT_TRUE(md_writer_ok(&fw.w));
    EXPECT_EQ((size_t)10, fw.passed);

    md_write_cstr(&fw.w, "0123456789abcdef");
    EXPECT_FALSE(md_writer_ok(&fw.w));

    // The error sticks: what follows is discarded, without reaching the sink, but is still counted
    md_write_cstr(&fw.w, "after");
    md_writef(&fw.w, "%s", "a format which does not fit");
    md_writer_flush(&fw.w);
    EXPECT_FALSE(md_writer_ok(&fw.w));
    EXPECT_EQ((size_t)10, fw.passed);
    EXPECT_EQ((size_t)(10 + 16 + 5 + 27), md_writer_tell(&fw.w));
}
