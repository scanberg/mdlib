#include "utest.h"
#include <core/md_allocator.h>
#include <core/md_base64.h>
#include <core/md_common.h>

#include <string.h>
#include <stdint.h>

UTEST(base64, test) {
    const char* test_strings[] = {
        "Hola, no habla español!",
        "This is a test string with some bogus...",
        "1",
        "FOUR"
    };

    char enc_buf[256];
    char dec_buf[256];

    for (int64_t i = 0; i < ARRAY_SIZE(test_strings); ++i) {
        int str_len = (int)strlen(test_strings[i]);
        int enc_len = md_base64_encode(enc_buf, test_strings[i], str_len);
        int dec_len = md_base64_decode(dec_buf, enc_buf, enc_len);
        EXPECT_EQ(str_len, dec_len);
        EXPECT_STRNEQ(test_strings[i], dec_buf, str_len);
    }
}