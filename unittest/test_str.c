#include "utest.h"

#include <core/md_str.h>
#include <core/md_os.h>
#include <core/md_allocator.h>
#include <float.h>

UTEST(str, extract_line) {
    str_t str = STR_INIT(
        "this is some text\n"
        "this is line 2\n"
        "\n"
        "\r\n"
        "}\n"
        "this is the end"
    );

    str_t line;
    
    EXPECT_TRUE(str_extract_line(&line, &str));
    EXPECT_STRNEQ("this is some text", line.ptr, line.len);
    
    EXPECT_TRUE(str_extract_line(&line, &str));
    EXPECT_STRNEQ("this is line 2", line.ptr, line.len);
    
    EXPECT_TRUE(str_extract_line(&line, &str));
    EXPECT_STRNEQ("", line.ptr, line.len);
    
    EXPECT_TRUE(str_extract_line(&line, &str));
    EXPECT_STRNEQ("", line.ptr, line.len);

    EXPECT_TRUE(str_extract_line(&line, &str));
    EXPECT_STRNEQ("}", line.ptr, line.len);

    EXPECT_TRUE(str_extract_line(&line, &str));
    EXPECT_STRNEQ("this is the end", line.ptr, line.len);

    EXPECT_FALSE(str_extract_line(&line, &str));
}

UTEST(str, edit_distance) {
    int dist;

    dist = str_edit_distance(STR_LIT("kitten"), STR_LIT("sitting"));
    EXPECT_EQ(3, dist);
    
    dist = str_edit_distance(STR_LIT("rosettacode"), STR_LIT("raisethysword"));
    EXPECT_EQ(8, dist);

    dist = str_edit_distance(STR_LIT(""), STR_LIT("something"));
    EXPECT_EQ(9, dist);
}

UTEST(str, count_equal_chars) {
    int count;
    count = str_count_equal_chars(STR_LIT("kitten"), STR_LIT("kittenz"));
    EXPECT_EQ(6, count);

    count = str_count_equal_chars(STR_LIT("kitten"), STR_LIT("sitting"));
    EXPECT_EQ(0, count);

    count = str_count_equal_chars(
        STR_LIT("/mnt/e/git/viamd/ext/mdlib/test_data/dir/subdir"),
        STR_LIT("/mnt/e/git/viamd/ext/mdlib/test_data/40-40-2-ddba-dyna.xmol"));
    EXPECT_EQ(sizeof("/mnt/e/git/viamd/ext/mdlib/test_data/") - 1, count);
}

UTEST(str, find_str) {
    size_t loc = SIZE_MAX;

    EXPECT_TRUE(str_find_str(&loc, STR_LIT("hello world"), STR_LIT("world")));
    EXPECT_EQ(loc, 6);
    EXPECT_TRUE(str_find_str(&loc, STR_LIT("hello"), STR_LIT("hello")));
    EXPECT_EQ(loc, 0);
    EXPECT_TRUE(str_find_str(&loc, STR_LIT("hello"), STR_LIT("o")));
    EXPECT_EQ(loc, 4);

    // The match starts inside an earlier partial match
    EXPECT_TRUE(str_find_str(&loc, STR_LIT("aab"), STR_LIT("ab")));
    EXPECT_EQ(loc, 1);
    EXPECT_TRUE(str_find_str(&loc, STR_LIT("ababac"), STR_LIT("abac")));
    EXPECT_EQ(loc, 2);
    EXPECT_TRUE(str_find_str(&loc, STR_LIT("<a id=\"x\"></a>"), STR_LIT("id=\"")));
    EXPECT_EQ(loc, 3);

    // First occurrence
    EXPECT_TRUE(str_find_str(&loc, STR_LIT("xyxy"), STR_LIT("xy")));
    EXPECT_EQ(loc, 0);

    // Nothing to find
    EXPECT_FALSE(str_find_str(&loc, STR_LIT("hello"), STR_LIT("world")));
    EXPECT_FALSE(str_find_str(&loc, STR_LIT("hell"),  STR_LIT("hello")));   // needle longer than haystack
    EXPECT_FALSE(str_find_str(&loc, STR_LIT("hello"), STR_LIT("")));
    EXPECT_FALSE(str_find_str(&loc, STR_LIT(""),      STR_LIT("a")));
    EXPECT_FALSE(str_find_str(&loc, STR_LIT("abc"),   STR_LIT("bcd")));     // would run past the end

    // loc is optional
    EXPECT_TRUE(str_find_str(NULL, STR_LIT("abc"), STR_LIT("bc")));
}

UTEST(str, copy_to_char_buf) {
    char buf[8];

    EXPECT_EQ(str_copy_to_char_buf(buf, sizeof(buf), STR_LIT("abc")), 3);
    EXPECT_STREQ("abc", buf);

    // Truncated to fit, and still zero terminated
    EXPECT_EQ(str_copy_to_char_buf(buf, sizeof(buf), STR_LIT("0123456789")), 7);
    EXPECT_STREQ("0123456", buf);

    // An empty string clears the buffer instead of leaving the previous contents
    EXPECT_EQ(str_copy_to_char_buf(buf, sizeof(buf), STR_LIT("")), 0);
    EXPECT_STREQ("", buf);
    str_copy_to_char_buf(buf, sizeof(buf), STR_LIT("abc"));
    EXPECT_EQ(str_copy_to_char_buf(buf, sizeof(buf), (str_t){0}), 0);
    EXPECT_STREQ("", buf);

    // A capacity of one only has room for the terminator, a capacity of zero writes nothing
    buf[0] = 'x';
    EXPECT_EQ(str_copy_to_char_buf(buf, 1, STR_LIT("abc")), 0);
    EXPECT_EQ(buf[0], '\0');
    buf[0] = 'x';
    EXPECT_EQ(str_copy_to_char_buf(buf, 0, STR_LIT("abc")), 0);
    EXPECT_EQ(buf[0], 'x');
}

// A str_t is a view: an empty one is the empty string wherever it points, as with std::string_view
UTEST(str, eq_empty) {
    const str_t null_str  = {0};
    const str_t lit       = STR_LIT("");
    const str_t slice     = str_substr(STR_LIT("abc"), 1, 0);   // empty, pointing into another string
    const str_t slice2    = str_substr(STR_LIT("xyz"), 2, 0);

    EXPECT_TRUE(str_eq(null_str, null_str));
    EXPECT_TRUE(str_eq(null_str, lit));
    EXPECT_TRUE(str_eq(lit, null_str));
    EXPECT_TRUE(str_eq(slice, slice2));
    EXPECT_TRUE(str_eq(slice, null_str));
    EXPECT_FALSE(str_eq(null_str, STR_LIT("a")));

    EXPECT_TRUE(str_eq_ignore_case(null_str, lit));
    EXPECT_TRUE(str_eq_n(null_str, lit, 4));
    EXPECT_TRUE(str_eq_n_ignore_case(slice, null_str, 4));

    EXPECT_TRUE(str_eq_cstr(null_str, ""));
    EXPECT_TRUE(str_eq_cstr(lit, ""));
    EXPECT_TRUE(str_eq_cstr(null_str, NULL));          // a NULL cstr is "", as for str_from_cstr
    EXPECT_FALSE(str_eq_cstr(STR_LIT("a"), NULL));
    EXPECT_FALSE(str_eq_cstr(null_str, "a"));
    EXPECT_TRUE(str_eq_cstr_ignore_case(slice, ""));
    EXPECT_TRUE(str_eq_cstr_n(null_str, "", 3));

    // A copy of an empty string is equal to it ({0} comes back from str_copy)
    md_allocator_i* alloc = md_get_heap_allocator();
    str_t copy = str_copy(lit, alloc);
    EXPECT_TRUE(str_eq(copy, lit));
    str_free(copy, alloc);
}

UTEST(str, eq) {
    EXPECT_TRUE(str_eq(STR_LIT("abc"), STR_LIT("abc")));
    EXPECT_FALSE(str_eq(STR_LIT("abc"), STR_LIT("abd")));
    EXPECT_FALSE(str_eq(STR_LIT("abc"), STR_LIT("ab")));
    EXPECT_TRUE(str_eq(str_substr(STR_LIT("xabcx"), 1, 3), STR_LIT("abc")));

    EXPECT_TRUE(str_eq_ignore_case(STR_LIT("HeLLo"), STR_LIT("hello")));
    EXPECT_FALSE(str_eq_ignore_case(STR_LIT("hello"), STR_LIT("hell")));

    EXPECT_TRUE(str_eq_cstr(STR_LIT("abc"), "abc"));
    EXPECT_FALSE(str_eq_cstr(STR_LIT("abc"), "abcd"));
    EXPECT_FALSE(str_eq_cstr(STR_LIT("abc"), "ab"));
    EXPECT_TRUE(str_eq_cstr(str_substr(STR_LIT("abcdef"), 0, 3), "abc"));   // not zero terminated
    EXPECT_TRUE(str_eq_cstr_ignore_case(STR_LIT("ATOM"), "atom"));
    EXPECT_FALSE(str_eq_cstr_ignore_case(STR_LIT("ATOM"), "atoms"));
}

// The _n variants compare the first n characters of both: with n = strlen(prefix) they test for a prefix,
// which a string shorter than the prefix does not have
UTEST(str, eq_n) {
    EXPECT_TRUE(str_eq_n(STR_LIT("ATOM      1"), STR_LIT("ATOM"), 4));
    EXPECT_TRUE(str_eq_n(STR_LIT("abcdef"), STR_LIT("abcxyz"), 3));
    EXPECT_FALSE(str_eq_n(STR_LIT("abcdef"), STR_LIT("abcxyz"), 4));
    EXPECT_FALSE(str_eq_n(STR_LIT("ab"), STR_LIT("abc"), 3));
    EXPECT_TRUE(str_eq_n(STR_LIT("ab"), STR_LIT("ab"), 3));       // both end within n, and are equal
    EXPECT_TRUE(str_eq_n(STR_LIT("abc"), STR_LIT("xyz"), 0));

    EXPECT_TRUE(str_eq_cstr_n(STR_LIT("ENDMDL"), "ENDMDL", 6));
    EXPECT_TRUE(str_eq_cstr_n(STR_LIT("ENDMDL      2"), "ENDMDL", 6));
    EXPECT_FALSE(str_eq_cstr_n(STR_LIT("END"), "ENDMDL", 6));
    EXPECT_FALSE(str_eq_cstr_n(STR_LIT("C"), "C3", 2));
    EXPECT_TRUE(str_eq_cstr_n(STR_LIT("C3'"), "C3", 2));
    EXPECT_TRUE(str_eq_cstr_n(STR_LIT("lj/cut"), "lj", 2));
    EXPECT_FALSE(str_eq_cstr_n(STR_LIT("REMARK"), "REMARK 350", 10));
    EXPECT_TRUE(str_eq_cstr_n(STR_LIT("ab"), "ab", 5));             // both end within n

    EXPECT_TRUE(str_eq_n_ignore_case(STR_LIT("Vibration 1"), STR_LIT("VIBRATION"), 9));
    EXPECT_TRUE(str_eq_cstr_n_ignore_case(STR_LIT("Time (ps)"), "time", 4));
    EXPECT_FALSE(str_eq_cstr_n_ignore_case(STR_LIT("Tim"), "time", 4));
}

UTEST(str, cmp_lex) {
    EXPECT_EQ(str_cmp_lex(STR_LIT("abc"), STR_LIT("abc")), 0);
    EXPECT_EQ(str_cmp_lex(STR_LIT("abc"), STR_LIT("abd")), -1);
    EXPECT_EQ(str_cmp_lex(STR_LIT("abd"), STR_LIT("abc")), 1);
    EXPECT_EQ(str_cmp_lex(STR_LIT("ab"),  STR_LIT("abc")), -1);    // a prefix first
    EXPECT_EQ(str_cmp_lex(STR_LIT("abc"), STR_LIT("ab")),  1);
    EXPECT_EQ(str_cmp_lex((str_t){0},     STR_LIT("")),    0);
    EXPECT_EQ(str_cmp_lex((str_t){0},     STR_LIT("a")),  -1);
    // Characters order as unsigned bytes: UTF-8 (lead byte >= 0x80) after ASCII, whatever the signedness of char
    EXPECT_EQ(str_cmp_lex(STR_LIT("z"), STR_LIT("\xC3\xA5")), -1);
}

UTEST(str, begins_ends_with) {
    EXPECT_TRUE(str_begins_with(STR_LIT("hello"), STR_LIT("he")));
    EXPECT_FALSE(str_begins_with(STR_LIT("he"), STR_LIT("hello")));
    EXPECT_TRUE(str_begins_with(STR_LIT("hello"), (str_t){0}));
    EXPECT_TRUE(str_begins_with((str_t){0}, STR_LIT("")));
    EXPECT_TRUE(str_ends_with(STR_LIT("file.xtc"), STR_LIT(".xtc")));
    EXPECT_FALSE(str_ends_with(STR_LIT("xtc"), STR_LIT(".xtc")));
    EXPECT_TRUE(str_ends_with(STR_LIT("file"), (str_t){0}));
}

UTEST(str, find_char) {
    size_t loc = SIZE_MAX;
    EXPECT_TRUE(str_find_char(&loc, STR_LIT("a=b=c"), '='));
    EXPECT_EQ(loc, 1);
    EXPECT_TRUE(str_rfind_char(&loc, STR_LIT("a=b=c"), '='));
    EXPECT_EQ(loc, 3);
    EXPECT_TRUE(str_rfind_char(&loc, STR_LIT("abc"), 'a'));
    EXPECT_EQ(loc, 0);
    EXPECT_FALSE(str_find_char(&loc, STR_LIT("abc"), 'x'));
    EXPECT_FALSE(str_rfind_char(&loc, STR_LIT("abc"), 'x'));
    EXPECT_FALSE(str_find_char(&loc, (str_t){0}, 'a'));
    EXPECT_FALSE(str_rfind_char(&loc, (str_t){0}, 'a'));
    EXPECT_TRUE(str_find_char(NULL, STR_LIT("abc"), 'c'));

    // A byte above 0x7F is found whether it is passed as a char or as an unsigned char value
    const str_t utf8 = STR_LIT("x\xC3\xA5y");
    EXPECT_TRUE(str_find_char(&loc, utf8, (char)0xC3));
    EXPECT_EQ(loc, 1);
    EXPECT_TRUE(str_find_char(&loc, utf8, 0xA5));
    EXPECT_EQ(loc, 2);
    EXPECT_TRUE(str_rfind_char(&loc, utf8, 0xC3));
    EXPECT_EQ(loc, 1);
}

UTEST(str, skip_peek_line) {
    str_t str = STR_LIT("first\r\nsecond\nthird");
    str_t line;

    EXPECT_TRUE(str_peek_line(&line, &str));
    EXPECT_TRUE(str_eq_cstr(line, "first"));
    EXPECT_TRUE(str_skip_line(&str));
    EXPECT_TRUE(str_extract_line(&line, &str));
    EXPECT_TRUE(str_eq_cstr(line, "second"));
    EXPECT_TRUE(str_peek_line(&line, &str));
    EXPECT_TRUE(str_eq_cstr(line, "third"));
    EXPECT_TRUE(str_skip_line(&str));
    EXPECT_FALSE(str_skip_line(&str));
    EXPECT_FALSE(str_peek_line(&line, &str));

    // A '\r' is only taken off together with a '\n'
    str = STR_LIT("end\r");
    EXPECT_TRUE(str_extract_line(&line, &str));
    EXPECT_TRUE(str_eq_cstr(line, "end\r"));
}

UTEST(str, swap) {
    str_t a = STR_LIT("a");
    str_t b = STR_LIT("bb");
    str_swap(&a, &b);
    EXPECT_TRUE(str_eq_cstr(a, "bb"));
    EXPECT_TRUE(str_eq_cstr(b, "a"));
}
