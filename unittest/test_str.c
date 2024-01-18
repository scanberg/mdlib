#include "utest.h"

#include <core/md_str.h>
#include <core/md_parse.h>
#include <core/md_os.h>
#include <float.h>

UTEST(str, parse_int) {
    str_t test_data[] = {
        STR_LIT("1023"),
        STR_LIT("-248"),
        STR_LIT("1232326745"),
        STR_LIT("1"),
        STR_LIT("0"),
        STR_LIT("-0"),
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
    const str_t str = STR_LIT("128326746123");
    const int64_t num_iter = 1000000;
    int64_t acc = 0;

    md_timestamp_t t0 = md_time_current();
    for (int64_t i = 0; i < num_iter; ++i) {
        acc += atol(str.ptr);
    }
    md_timestamp_t t1 = md_time_current();    
    for (int64_t i = 0; i < num_iter; ++i) {
        acc += parse_int(str);
    }
    md_timestamp_t t2 = md_time_current();

    double t_atoi  = md_time_as_milliseconds(t1 - t0);
    double t_parse = md_time_as_milliseconds(t2 - t1);

    printf("Time to parse %iM int. atoi: %.3f ms, parse_int: %.3f ms, speedup: %.2f, %i\n", (int)(num_iter / 1000000), t_atoi, t_parse, t_atoi / t_parse, (int)acc);
}


UTEST(str, parse_float) {
    str_t test_data[] = {
        STR_LIT("1023.22311283798172389718923789172389"),
        STR_LIT("-248.273"),
        STR_LIT("0000000000000.273"),
        STR_LIT("1232326745e10"),
        STR_LIT("1.0e-29"),
        STR_LIT("0.02e+10"),
        STR_LIT("-0"),
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
    const str_t str = STR_LIT("-248.271233");
    const int64_t num_iter = 1000000;
    double acc = 0;

    md_timestamp_t t0 = md_time_current();
    for (int64_t i = 0; i < num_iter; ++i) {
        acc += atof(str.ptr);
    }
    md_timestamp_t t1 = md_time_current();    
    for (int64_t i = 0; i < num_iter; ++i) {
        acc += parse_float(str);
    }
    md_timestamp_t t2 = md_time_current();

    double t_atof  = md_time_as_milliseconds(t1 - t0);
    double t_parse = md_time_as_milliseconds(t2 - t1);

    printf("Time to parse %iM floats. atof: %.3f ms, parse_float: %.3f ms, speedup: %.2f, acc: %.1f\n", (int)(num_iter / 1000000), t_atof, t_parse, t_atof / t_parse, acc);
}

UTEST(str, extract_line) {
    str_t str = STR_LIT(
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