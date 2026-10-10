#include "utest.h"

#include <core/md_str.h>
#include <core/md_os.h>
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
