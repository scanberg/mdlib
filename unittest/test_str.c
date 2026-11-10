#include "utest.h"

#include <core/md_str.h>
#include <core/md_os.h>
#include <core/md_allocator.h>
#include <core/md_writer.h>
#include <md_types.h>
#include <float.h>

#include "comma_locale.h"

// What mdlib formats reads back with md_parse, which does not follow the locale either: under a decimal
// comma locale 1.5 is still written "1.5", by md_snprintf and by everything formatting through it
UTEST(str, printf_ignores_the_locale) {
    comma_locale_t loc;
    if (!comma_locale_enter(&loc)) {
        return;
    }

    char buf[128];
    // The C library follows the locale, or this would check nothing
    snprintf(buf, sizeof(buf), "%.2f", 1.5);
    EXPECT_STREQ("1,50", buf);

    EXPECT_EQ(4, md_snprintf(buf, sizeof(buf), "%.2f", 1.5));
    EXPECT_STREQ("1.50", buf);
    md_snprintf(buf, sizeof(buf), "%.9g|%.17g|%e|%g", (double)0.1f, 0.25, 1.5e-7, -2.5);
    EXPECT_STREQ("0.100000001|0.25|1.500000e-07|-2.5", buf);

    // Measured, and cut short, the way C99 has it
    EXPECT_EQ(4, md_snprintf(NULL, 0, "%.2f", 1.5));
    EXPECT_EQ(4, md_snprintf(buf, 3, "%.2f", 1.5));
    EXPECT_STREQ("1.", buf);

    md_allocator_i* alloc = md_get_heap_allocator();
    md_str_t str = md_str_printf(alloc, "%.3f", 2.25);
    EXPECT_TRUE(md_str_eq_cstr(str, "2.250"));
    md_str_free(str, alloc);

    md_str_writer_t sw = md_str_writer_create(alloc);
    md_str_writer_fmt(&sw, "%g,%g", 0.5, 1.5);
    EXPECT_TRUE(md_str_eq_cstr(md_str_writer_to_str(&sw), "0.5,1.5"));
    md_str_writer_free(&sw);

    comma_locale_leave(&loc);
}

UTEST(str, extract_line) {
    md_str_t str = MD_STR_INIT(
        "this is some text\n"
        "this is line 2\n"
        "\n"
        "\r\n"
        "}\n"
        "this is the end"
    );

    md_str_t line;
    
    EXPECT_TRUE(md_str_extract_line(&line, &str));
    EXPECT_STRNEQ("this is some text", line.ptr, line.len);
    
    EXPECT_TRUE(md_str_extract_line(&line, &str));
    EXPECT_STRNEQ("this is line 2", line.ptr, line.len);
    
    EXPECT_TRUE(md_str_extract_line(&line, &str));
    EXPECT_STRNEQ("", line.ptr, line.len);
    
    EXPECT_TRUE(md_str_extract_line(&line, &str));
    EXPECT_STRNEQ("", line.ptr, line.len);

    EXPECT_TRUE(md_str_extract_line(&line, &str));
    EXPECT_STRNEQ("}", line.ptr, line.len);

    EXPECT_TRUE(md_str_extract_line(&line, &str));
    EXPECT_STRNEQ("this is the end", line.ptr, line.len);

    EXPECT_FALSE(md_str_extract_line(&line, &str));
}

UTEST(str, edit_distance) {
    int dist;

    dist = md_str_edit_distance(MD_STR_LIT("kitten"), MD_STR_LIT("sitting"));
    EXPECT_EQ(3, dist);
    
    dist = md_str_edit_distance(MD_STR_LIT("rosettacode"), MD_STR_LIT("raisethysword"));
    EXPECT_EQ(8, dist);

    dist = md_str_edit_distance(MD_STR_LIT(""), MD_STR_LIT("something"));
    EXPECT_EQ(9, dist);
}

UTEST(str, count_equal_chars) {
    int count;
    count = md_str_count_equal_chars(MD_STR_LIT("kitten"), MD_STR_LIT("kittenz"));
    EXPECT_EQ(6, count);

    count = md_str_count_equal_chars(MD_STR_LIT("kitten"), MD_STR_LIT("sitting"));
    EXPECT_EQ(0, count);

    count = md_str_count_equal_chars(
        MD_STR_LIT("/mnt/e/git/viamd/ext/mdlib/test_data/dir/subdir"),
        MD_STR_LIT("/mnt/e/git/viamd/ext/mdlib/test_data/40-40-2-ddba-dyna.xmol"));
    EXPECT_EQ(sizeof("/mnt/e/git/viamd/ext/mdlib/test_data/") - 1, count);
}

UTEST(str, find_str) {
    size_t loc = SIZE_MAX;

    EXPECT_TRUE(md_str_find_str(&loc, MD_STR_LIT("hello world"), MD_STR_LIT("world")));
    EXPECT_EQ(loc, 6);
    EXPECT_TRUE(md_str_find_str(&loc, MD_STR_LIT("hello"), MD_STR_LIT("hello")));
    EXPECT_EQ(loc, 0);
    EXPECT_TRUE(md_str_find_str(&loc, MD_STR_LIT("hello"), MD_STR_LIT("o")));
    EXPECT_EQ(loc, 4);

    // The match starts inside an earlier partial match
    EXPECT_TRUE(md_str_find_str(&loc, MD_STR_LIT("aab"), MD_STR_LIT("ab")));
    EXPECT_EQ(loc, 1);
    EXPECT_TRUE(md_str_find_str(&loc, MD_STR_LIT("ababac"), MD_STR_LIT("abac")));
    EXPECT_EQ(loc, 2);
    EXPECT_TRUE(md_str_find_str(&loc, MD_STR_LIT("<a id=\"x\"></a>"), MD_STR_LIT("id=\"")));
    EXPECT_EQ(loc, 3);

    // First occurrence
    EXPECT_TRUE(md_str_find_str(&loc, MD_STR_LIT("xyxy"), MD_STR_LIT("xy")));
    EXPECT_EQ(loc, 0);

    // Nothing to find
    EXPECT_FALSE(md_str_find_str(&loc, MD_STR_LIT("hello"), MD_STR_LIT("world")));
    EXPECT_FALSE(md_str_find_str(&loc, MD_STR_LIT("hell"),  MD_STR_LIT("hello")));   // needle longer than haystack
    EXPECT_FALSE(md_str_find_str(&loc, MD_STR_LIT("hello"), MD_STR_LIT("")));
    EXPECT_FALSE(md_str_find_str(&loc, MD_STR_LIT(""),      MD_STR_LIT("a")));
    EXPECT_FALSE(md_str_find_str(&loc, MD_STR_LIT("abc"),   MD_STR_LIT("bcd")));     // would run past the end

    // loc is optional
    EXPECT_TRUE(md_str_find_str(NULL, MD_STR_LIT("abc"), MD_STR_LIT("bc")));
}

UTEST(str, copy_to_char_buf) {
    char buf[8];

    EXPECT_EQ(md_str_copy_to_char_buf(buf, sizeof(buf), MD_STR_LIT("abc")), 3);
    EXPECT_STREQ("abc", buf);

    // Truncated to fit, and still zero terminated
    EXPECT_EQ(md_str_copy_to_char_buf(buf, sizeof(buf), MD_STR_LIT("0123456789")), 7);
    EXPECT_STREQ("0123456", buf);

    // An empty string clears the buffer instead of leaving the previous contents
    EXPECT_EQ(md_str_copy_to_char_buf(buf, sizeof(buf), MD_STR_LIT("")), 0);
    EXPECT_STREQ("", buf);
    md_str_copy_to_char_buf(buf, sizeof(buf), MD_STR_LIT("abc"));
    EXPECT_EQ(md_str_copy_to_char_buf(buf, sizeof(buf), (md_str_t){0}), 0);
    EXPECT_STREQ("", buf);

    // A capacity of one only has room for the terminator, a capacity of zero writes nothing
    buf[0] = 'x';
    EXPECT_EQ(md_str_copy_to_char_buf(buf, 1, MD_STR_LIT("abc")), 0);
    EXPECT_EQ(buf[0], '\0');
    buf[0] = 'x';
    EXPECT_EQ(md_str_copy_to_char_buf(buf, 0, MD_STR_LIT("abc")), 0);
    EXPECT_EQ(buf[0], 'x');
}

// A md_str_t is a view: an empty one is the empty string wherever it points, as with std::string_view
UTEST(str, eq_empty) {
    const md_str_t null_str  = {0};
    const md_str_t lit       = MD_STR_LIT("");
    const md_str_t slice     = md_str_substr(MD_STR_LIT("abc"), 1, 0);   // empty, pointing into another string
    const md_str_t slice2    = md_str_substr(MD_STR_LIT("xyz"), 2, 0);

    EXPECT_TRUE(md_str_eq(null_str, null_str));
    EXPECT_TRUE(md_str_eq(null_str, lit));
    EXPECT_TRUE(md_str_eq(lit, null_str));
    EXPECT_TRUE(md_str_eq(slice, slice2));
    EXPECT_TRUE(md_str_eq(slice, null_str));
    EXPECT_FALSE(md_str_eq(null_str, MD_STR_LIT("a")));

    EXPECT_TRUE(md_str_eq_ignore_case(null_str, lit));
    EXPECT_TRUE(md_str_eq_n(null_str, lit, 4));
    EXPECT_TRUE(md_str_eq_n_ignore_case(slice, null_str, 4));

    EXPECT_TRUE(md_str_eq_cstr(null_str, ""));
    EXPECT_TRUE(md_str_eq_cstr(lit, ""));
    EXPECT_TRUE(md_str_eq_cstr(null_str, NULL));          // a NULL cstr is "", as for md_str_from_cstr
    EXPECT_FALSE(md_str_eq_cstr(MD_STR_LIT("a"), NULL));
    EXPECT_FALSE(md_str_eq_cstr(null_str, "a"));
    EXPECT_TRUE(md_str_eq_cstr_ignore_case(slice, ""));
    EXPECT_TRUE(md_str_eq_cstr_n(null_str, "", 3));

    // A copy of an empty string is equal to it ({0} comes back from md_str_copy)
    md_allocator_i* alloc = md_get_heap_allocator();
    md_str_t copy = md_str_copy(lit, alloc);
    EXPECT_TRUE(md_str_eq(copy, lit));
    md_str_free(copy, alloc);
}

UTEST(str, eq) {
    EXPECT_TRUE(md_str_eq(MD_STR_LIT("abc"), MD_STR_LIT("abc")));
    EXPECT_FALSE(md_str_eq(MD_STR_LIT("abc"), MD_STR_LIT("abd")));
    EXPECT_FALSE(md_str_eq(MD_STR_LIT("abc"), MD_STR_LIT("ab")));
    EXPECT_TRUE(md_str_eq(md_str_substr(MD_STR_LIT("xabcx"), 1, 3), MD_STR_LIT("abc")));

    EXPECT_TRUE(md_str_eq_ignore_case(MD_STR_LIT("HeLLo"), MD_STR_LIT("hello")));
    EXPECT_FALSE(md_str_eq_ignore_case(MD_STR_LIT("hello"), MD_STR_LIT("hell")));

    EXPECT_TRUE(md_str_eq_cstr(MD_STR_LIT("abc"), "abc"));
    EXPECT_FALSE(md_str_eq_cstr(MD_STR_LIT("abc"), "abcd"));
    EXPECT_FALSE(md_str_eq_cstr(MD_STR_LIT("abc"), "ab"));
    EXPECT_TRUE(md_str_eq_cstr(md_str_substr(MD_STR_LIT("abcdef"), 0, 3), "abc"));   // not zero terminated
    EXPECT_TRUE(md_str_eq_cstr_ignore_case(MD_STR_LIT("ATOM"), "atom"));
    EXPECT_FALSE(md_str_eq_cstr_ignore_case(MD_STR_LIT("ATOM"), "atoms"));
}

// The _n variants compare the first n characters of both: with n = strlen(prefix) they test for a prefix,
// which a string shorter than the prefix does not have
UTEST(str, eq_n) {
    EXPECT_TRUE(md_str_eq_n(MD_STR_LIT("ATOM      1"), MD_STR_LIT("ATOM"), 4));
    EXPECT_TRUE(md_str_eq_n(MD_STR_LIT("abcdef"), MD_STR_LIT("abcxyz"), 3));
    EXPECT_FALSE(md_str_eq_n(MD_STR_LIT("abcdef"), MD_STR_LIT("abcxyz"), 4));
    EXPECT_FALSE(md_str_eq_n(MD_STR_LIT("ab"), MD_STR_LIT("abc"), 3));
    EXPECT_TRUE(md_str_eq_n(MD_STR_LIT("ab"), MD_STR_LIT("ab"), 3));       // both end within n, and are equal
    EXPECT_TRUE(md_str_eq_n(MD_STR_LIT("abc"), MD_STR_LIT("xyz"), 0));

    EXPECT_TRUE(md_str_eq_cstr_n(MD_STR_LIT("ENDMDL"), "ENDMDL", 6));
    EXPECT_TRUE(md_str_eq_cstr_n(MD_STR_LIT("ENDMDL      2"), "ENDMDL", 6));
    EXPECT_FALSE(md_str_eq_cstr_n(MD_STR_LIT("END"), "ENDMDL", 6));
    EXPECT_FALSE(md_str_eq_cstr_n(MD_STR_LIT("C"), "C3", 2));
    EXPECT_TRUE(md_str_eq_cstr_n(MD_STR_LIT("C3'"), "C3", 2));
    EXPECT_TRUE(md_str_eq_cstr_n(MD_STR_LIT("lj/cut"), "lj", 2));
    EXPECT_FALSE(md_str_eq_cstr_n(MD_STR_LIT("REMARK"), "REMARK 350", 10));
    EXPECT_TRUE(md_str_eq_cstr_n(MD_STR_LIT("ab"), "ab", 5));             // both end within n

    EXPECT_TRUE(md_str_eq_n_ignore_case(MD_STR_LIT("Vibration 1"), MD_STR_LIT("VIBRATION"), 9));
    EXPECT_TRUE(md_str_eq_cstr_n_ignore_case(MD_STR_LIT("Time (ps)"), "time", 4));
    EXPECT_FALSE(md_str_eq_cstr_n_ignore_case(MD_STR_LIT("Tim"), "time", 4));
}

UTEST(str, cmp_lex) {
    EXPECT_EQ(md_str_cmp_lex(MD_STR_LIT("abc"), MD_STR_LIT("abc")), 0);
    EXPECT_EQ(md_str_cmp_lex(MD_STR_LIT("abc"), MD_STR_LIT("abd")), -1);
    EXPECT_EQ(md_str_cmp_lex(MD_STR_LIT("abd"), MD_STR_LIT("abc")), 1);
    EXPECT_EQ(md_str_cmp_lex(MD_STR_LIT("ab"),  MD_STR_LIT("abc")), -1);    // a prefix first
    EXPECT_EQ(md_str_cmp_lex(MD_STR_LIT("abc"), MD_STR_LIT("ab")),  1);
    EXPECT_EQ(md_str_cmp_lex((md_str_t){0},     MD_STR_LIT("")),    0);
    EXPECT_EQ(md_str_cmp_lex((md_str_t){0},     MD_STR_LIT("a")),  -1);
    // Characters order as unsigned bytes: UTF-8 (lead byte >= 0x80) after ASCII, whatever the signedness of char
    EXPECT_EQ(md_str_cmp_lex(MD_STR_LIT("z"), MD_STR_LIT("\xC3\xA5")), -1);
}

UTEST(str, begins_ends_with) {
    EXPECT_TRUE(md_str_begins_with(MD_STR_LIT("hello"), MD_STR_LIT("he")));
    EXPECT_FALSE(md_str_begins_with(MD_STR_LIT("he"), MD_STR_LIT("hello")));
    EXPECT_TRUE(md_str_begins_with(MD_STR_LIT("hello"), (md_str_t){0}));
    EXPECT_TRUE(md_str_begins_with((md_str_t){0}, MD_STR_LIT("")));
    EXPECT_TRUE(md_str_ends_with(MD_STR_LIT("file.xtc"), MD_STR_LIT(".xtc")));
    EXPECT_FALSE(md_str_ends_with(MD_STR_LIT("xtc"), MD_STR_LIT(".xtc")));
    EXPECT_TRUE(md_str_ends_with(MD_STR_LIT("file"), (md_str_t){0}));
}

UTEST(str, find_char) {
    size_t loc = SIZE_MAX;
    EXPECT_TRUE(md_str_find_char(&loc, MD_STR_LIT("a=b=c"), '='));
    EXPECT_EQ(loc, 1);
    EXPECT_TRUE(md_str_rfind_char(&loc, MD_STR_LIT("a=b=c"), '='));
    EXPECT_EQ(loc, 3);
    EXPECT_TRUE(md_str_rfind_char(&loc, MD_STR_LIT("abc"), 'a'));
    EXPECT_EQ(loc, 0);
    EXPECT_FALSE(md_str_find_char(&loc, MD_STR_LIT("abc"), 'x'));
    EXPECT_FALSE(md_str_rfind_char(&loc, MD_STR_LIT("abc"), 'x'));
    EXPECT_FALSE(md_str_find_char(&loc, (md_str_t){0}, 'a'));
    EXPECT_FALSE(md_str_rfind_char(&loc, (md_str_t){0}, 'a'));
    EXPECT_TRUE(md_str_find_char(NULL, MD_STR_LIT("abc"), 'c'));

    // A byte above 0x7F is found whether it is passed as a char or as an unsigned char value
    const md_str_t utf8 = MD_STR_LIT("x\xC3\xA5y");
    EXPECT_TRUE(md_str_find_char(&loc, utf8, (char)0xC3));
    EXPECT_EQ(loc, 1);
    EXPECT_TRUE(md_str_find_char(&loc, utf8, 0xA5));
    EXPECT_EQ(loc, 2);
    EXPECT_TRUE(md_str_rfind_char(&loc, utf8, 0xC3));
    EXPECT_EQ(loc, 1);
}

UTEST(str, skip_peek_line) {
    md_str_t str = MD_STR_LIT("first\r\nsecond\nthird");
    md_str_t line;

    EXPECT_TRUE(md_str_peek_line(&line, &str));
    EXPECT_TRUE(md_str_eq_cstr(line, "first"));
    EXPECT_TRUE(md_str_skip_line(&str));
    EXPECT_TRUE(md_str_extract_line(&line, &str));
    EXPECT_TRUE(md_str_eq_cstr(line, "second"));
    EXPECT_TRUE(md_str_peek_line(&line, &str));
    EXPECT_TRUE(md_str_eq_cstr(line, "third"));
    EXPECT_TRUE(md_str_skip_line(&str));
    EXPECT_FALSE(md_str_skip_line(&str));
    EXPECT_FALSE(md_str_peek_line(&line, &str));

    // A '\r' is only taken off together with a '\n'
    str = MD_STR_LIT("end\r");
    EXPECT_TRUE(md_str_extract_line(&line, &str));
    EXPECT_TRUE(md_str_eq_cstr(line, "end\r"));
}

UTEST(str, swap) {
    md_str_t a = MD_STR_LIT("a");
    md_str_t b = MD_STR_LIT("bb");
    md_str_swap(&a, &b);
    EXPECT_TRUE(md_str_eq_cstr(a, "bb"));
    EXPECT_TRUE(md_str_eq_cstr(b, "a"));
}

// A label holds up to 8 characters (md_label16_t 16), zero after them: its length is where the first zero byte is,
// all of it when there is none, and two labels are equal when their bytes are
UTEST(str, labels) {
    static const char text[] = "ABCDEFGHIJKLMNOPQRS";
    for (size_t n = 0; n <= 10; ++n) {
        const md_label_t l = make_label((md_str_t){text, n});
        EXPECT_EQ(MD_MIN(n, 8), md_label_len(&l));
        EXPECT_TRUE(md_str_eq(md_label_str(&l), (md_str_t){text, MD_MIN(n, 8)}));
        EXPECT_EQ(n == 0, label_empty(l));
    }
    for (size_t n = 0; n <= 18; ++n) {
        const md_label16_t l = make_label16((md_str_t){text, n});
        EXPECT_EQ(MD_MIN(n, 16), md_label16_len(&l));
        EXPECT_TRUE(md_str_eq(md_label16_str(&l), (md_str_t){text, MD_MIN(n, 16)}));
    }
    // Bytes past ASCII count as characters
    const md_label_t u = make_label(MD_STR_LIT("C\xc3\xa5\x80\xff"));
    EXPECT_EQ(5, md_label_len(&u));

    EXPECT_TRUE(md_label_eq(make_label(MD_STR_LIT("CA")), make_label(MD_STR_LIT("CA"))));
    EXPECT_FALSE(md_label_eq(make_label(MD_STR_LIT("CA")), make_label(MD_STR_LIT("CA1"))));
    EXPECT_TRUE(md_label_eq(make_label(MD_STR_LIT("ABCDEFGH")), make_label(MD_STR_LIT("ABCDEFGHI"))));
    EXPECT_NE(md_label_key(make_label(MD_STR_LIT("CA"))), md_label_key(make_label(MD_STR_LIT("CB"))));
}

UTEST(str, extract_ext_of_the_file_name) {
    md_str_t ext = {0};
    EXPECT_TRUE(md_extract_ext(&ext, MD_STR_LIT("/data/run.1/traj.nc")));
    EXPECT_TRUE(md_str_eq_cstr(ext, "nc"));
    // A dot in a folder is not an extension of a file without one (Amber's mdout, prmtop)
    EXPECT_FALSE(md_extract_ext(&ext, MD_STR_LIT("/data/run.1/mdout")));
    EXPECT_FALSE(md_extract_ext(&ext, MD_STR_LIT("C:\\runs\\v1.2\\prmtop")));
    EXPECT_TRUE(md_extract_ext(&ext, MD_STR_LIT("C:\\runs\\v1.2\\sys.parm7")));
    EXPECT_TRUE(md_str_eq_cstr(ext, "parm7"));
}
