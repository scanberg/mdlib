#include "utest.h"

#include <core/md_str.h>
#include <core/md_os.h>
#include <float.h>

UTEST(str, parse_int) {
    str_t test_data[] = {
        STR("1023"),
        STR("-248"),
        STR("1232326745"),
        STR("1"),
        STR("0"),
        STR("-0"),
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
    const str_t str = STR("128326746123");
    const int64_t num_iter = 1000000;

    md_timestamp_t t0 = md_time_current();
    for (int64_t i = 0; i < 100000; ++i) {
        int64_t x = atol(str.ptr);
    }
    md_timestamp_t t1 = md_time_current();    
    for (int64_t i = 0; i < 100000; ++i) {
        int64_t x = parse_int(str);
    }
    md_timestamp_t t2 = md_time_current();

    double t_atoi  = md_time_as_milliseconds(t1 - t0);
    double t_parse = md_time_as_milliseconds(t2 - t1);

    printf("Time to parse %iM int. atoi: %.3f ms, parse_int: %.3f ms, speedup: %.2f\n", (int)(num_iter / 1000000), t_atoi, t_parse, t_atoi / t_parse);
}


UTEST(str, parse_float) {
    str_t test_data[] = {
        STR("1023.22311283798172389718923789172389"),
        STR("-248.273"),
        STR("0000000000000.273"),
        STR("1232326745e10"),
        STR("1.0e-29"),
        STR("0.02e+10"),
        STR("-0"),
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
    const str_t str = STR("-248.271233e-2");
    const int64_t num_iter = 1000000;

    md_timestamp_t t0 = md_time_current();
    for (int64_t i = 0; i < 100000; ++i) {
        double x = atof(str.ptr);
    }
    md_timestamp_t t1 = md_time_current();    
    for (int64_t i = 0; i < 100000; ++i) {
        double x = parse_float(str);
    }
    md_timestamp_t t2 = md_time_current();

    double t_atof  = md_time_as_milliseconds(t1 - t0);
    double t_parse = md_time_as_milliseconds(t2 - t1);

    printf("Time to parse %iM floats. atof: %.3f ms, parse_float: %.3f ms, speedup: %.2f\n", (int)(num_iter / 1000000), t_atof, t_parse, t_atof / t_parse);
}

UTEST(str, extract_line) {
    const str_t str = STR(
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
