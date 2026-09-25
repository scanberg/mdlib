#include "utest.h"

#include <core/md_str.h>
#include <core/md_parse.h>
#include <core/md_os.h>
#include <float.h>

UTEST(str, parse_int) {
    str_t test_data[] = {
        STR_INIT("1023"),
        STR_INIT("-248"),
        STR_INIT("1232326745"),
        STR_INIT("1"),
        STR_INIT("0"),
        STR_INIT("-0"),
    };

    int64_t ref_data[] = {
        1023, -248, 1232326745, 1, 0, 0,
    };

    for (int64_t i = 0; i < ARRAY_SIZE(test_data); ++i) {
        int64_t val = parse_int(test_data[i]);
        int64_t ref = ref_data[i];
        EXPECT_EQ(ref, val);
    }
}

UTEST(str, perf_int) {
    const str_t str = STR_INIT("128326746123");
    const int64_t num_iter = 1000000;
    int64_t acc = 0;

    md_tick_t t0 = md_tick_now();
    for (int64_t i = 0; i < num_iter; ++i) {
        acc += atol(str.ptr);
    }
    md_tick_t t1 = md_tick_now();    
    for (int64_t i = 0; i < num_iter; ++i) {
        acc += parse_int(str);
    }
    md_tick_t t2 = md_tick_now();

    double t_atoi  = md_tick_to_milliseconds(t1 - t0);
    double t_parse = md_tick_to_milliseconds(t2 - t1);

    printf("Time to parse %iM int. atoi: %.3f ms, parse_int: %.3f ms, speedup: %.2f, %i\n", (int)(num_iter / 1000000), t_atoi, t_parse, t_atoi / t_parse, (int)acc);
}


UTEST(str, parse_float) {
    str_t test_data[] = {
        STR_INIT("1023.22311283798172389718923789172389"),
        STR_INIT("-248.273"),
        STR_INIT("0000000000000.273"),
        STR_INIT("1232326745e10"),
        STR_INIT("1.0e-29"),
        STR_INIT("0.02e+10"),
        STR_INIT("-0"),
    };

    double ref_data[] = {
        1023.22311283798172389718923789172389,
        -248.273,
        0000000000000.273,
        1232326745e10,
        1.0e-16,
        0.02e+10,
        -0,
    };

    for (int64_t i = 0; i < ARRAY_SIZE(test_data); ++i) {
        double val = parse_float(test_data[i]);
        double ref = ref_data[i];
        EXPECT_NEAR(ref, val, 1e-12);
    }
}

UTEST(str, perf_float) {
    const str_t str = STR_INIT("-248.271233");
    const int64_t num_iter = 1000000;
    double acc = 0;

    md_tick_t t0 = md_tick_now();
    for (int64_t i = 0; i < num_iter; ++i) {
        acc += atof(str.ptr);
    }
    md_tick_t t1 = md_tick_now();    
    for (int64_t i = 0; i < num_iter; ++i) {
        acc += parse_float(str);
    }
    md_tick_t t2 = md_tick_now();

    double t_atof  = md_tick_to_milliseconds(t1 - t0);
    double t_parse = md_tick_to_milliseconds(t2 - t1);

    printf("Time to parse %iM floats. atof: %.3f ms, parse_float: %.3f ms, speedup: %.2f, acc: %.1f\n", (int)(num_iter / 1000000), t_atof, t_parse, t_atof / t_parse, acc);
}

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
