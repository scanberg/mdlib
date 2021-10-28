#include "utest.h"
#include <core/md_allocator.h>
#include <core/md_base64.h>
#include <core/md_common.h>

#include <string.h>
#include <stdint.h>

UTEST(base64, encode_decode) {
    const char* input[] = {
        "ñ",
        "This is a test string with some bogus...",
        "1",
        "FOUR"
    };

    const char* ref_enc[] = {
        "w7E=",
        "VGhpcyBpcyBhIHRlc3Qgc3RyaW5nIHdpdGggc29tZSBib2d1cy4uLg==",
        "MQ==",
        "Rk9VUg=="
    };

    char enc_buf[256];
    char dec_buf[256];

    for (int64_t i = 0; i < ARRAY_SIZE(input); ++i) {
        int str_len = (int)strlen(input[i]);
        int enc_len = md_base64_encode(enc_buf, input[i], str_len);

        EXPECT_EQ(enc_len, strlen(ref_enc[i]));
        EXPECT_STRNEQ(ref_enc[i], enc_buf, enc_len);

        int dec_len = md_base64_decode(dec_buf, enc_buf, enc_len);
        EXPECT_EQ(str_len, dec_len);
        EXPECT_STRNEQ(input[i], dec_buf, str_len);
    }
}

UTEST(base64, enc_with_new_line) {
    const char input[] = "This is a test string with some bogus...";
    const char ref_enc[] = "VGhpcyBpcyBhIHRlc3Qgc3RyaW5nIHdpdGggc29tZSBib2d1cy4uLg==";

    char enc_buf[256];
    char dec_buf[256];

    int input_len = sizeof(input)-1;
    int enc_len = md_base64_encode(enc_buf, input, input_len);

    EXPECT_EQ(enc_len, strlen(ref_enc));
    EXPECT_STRNEQ(ref_enc, enc_buf, enc_len);

    char enc_with_newline[256];
    int enc_with_newline_len = snprintf(enc_with_newline, sizeof(enc_with_newline), "%.*s\r\n%.*s", 30, enc_buf, enc_len - 30, enc_buf + 30);

    int dec_len = md_base64_decode(dec_buf, enc_with_newline, enc_with_newline_len);
    EXPECT_EQ(input_len, dec_len);
    EXPECT_STRNEQ(input, dec_buf, dec_len);
}

UTEST(base64, not_multiple_of_four) {
    const char enc[] = "VGhpyBp";
    char dec_buf[256];
    int dec_len = md_base64_decode(dec_buf, enc, (int)strlen(enc));
    EXPECT_EQ(0, dec_len);
}

UTEST(base64, invalid_character) {
    const char enc[] = "VGhpcyBp!ÅÄÖ";
    char dec_buf[256];
    int dec_len = md_base64_decode(dec_buf, enc, (int)strlen(enc));
    EXPECT_EQ(0, dec_len);
}