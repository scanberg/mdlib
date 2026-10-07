#include "ubench.h"

#include <core/md_str.h>
#include <core/md_parse.h>
#include <core/md_os.h>
#include <core/md_allocator.h>
#include <core/md_log.h>

#include <inttypes.h>
#include <stdlib.h>
#include <string.h>

UBENCH_EX(str, buffered_reader) {
    str_t path = STR_INIT(MD_BENCHMARK_DATA_DIR "/centered.gro");
    md_file_t file = {0};
    if (!md_file_open(&file, path, MD_FILE_READ)) {
        MD_LOG_ERROR("Could not open file '%.*s'", path.len, path.ptr);
        return;
    }
    const int64_t cap = MEGABYTES(1);
    char* buf = md_alloc(md_get_heap_allocator(), cap);
    
    md_buffered_reader_t reader = md_buffered_reader_from_file(buf, cap, file);

    UBENCH_SET_BYTES(md_file_size(file));

    UBENCH_DO_BENCHMARK() {
        md_file_seek(file, 0, MD_FILE_BEG);
        str_t line;
        while (md_buffered_reader_extract_line(&line, &reader)) {
            // do nothing
        }
    }

    md_free(md_get_heap_allocator(), buf, cap);
    md_file_close(&file);
}

UBENCH_EX(str, parse_int) {
    str_t str[] = {
        STR_INIT("1928123123123"),
        STR_INIT("1123    "),
        STR_INIT("19228123"),
        STR_INIT("1921238123"),
    };

    int64_t num_bytes = 0;
    for (int i = 0; i < (int)ARRAY_SIZE(str); ++i) {
        num_bytes += str[i].len;
    }
    UBENCH_SET_BYTES(num_bytes);

    size_t acc = 0;
    UBENCH_DO_BENCHMARK() {
		acc += parse_int(str[0]);
        acc += parse_int(str[1]);
        acc += parse_int(str[2]);
        acc += parse_int(str[3]);
    }
    UBENCH_DO_NOTHING(&acc);
}

// The C library on the same input, for reference
UBENCH_EX(str, strtoll) {
    const char* str[] = {
        "1928123123123",
        "1123    ",
        "19228123",
        "1921238123",
    };

    int64_t num_bytes = 0;
    for (int i = 0; i < (int)ARRAY_SIZE(str); ++i) {
        num_bytes += strlen(str[i]);
    }
    UBENCH_SET_BYTES(num_bytes);

    long long acc = 0;
    UBENCH_DO_BENCHMARK() {
        acc += strtoll(str[0], NULL, 10);
        acc += strtoll(str[1], NULL, 10);
        acc += strtoll(str[2], NULL, 10);
        acc += strtoll(str[3], NULL, 10);
    }
    UBENCH_DO_NOTHING(&acc);
}

UBENCH_EX(str, parse_float) {
    str_t str[] = {
        STR_INIT("1928123.2767"),
        STR_INIT("19.2    "),
        STR_INIT("12323   "),
        STR_INIT("0.000000"),
    };

    int64_t num_bytes = 0;
    for (int i = 0; i < (int)ARRAY_SIZE(str); ++i) {
        num_bytes += str[i].len;
    }
    UBENCH_SET_BYTES(num_bytes);

    double acc = 0;
    UBENCH_DO_BENCHMARK() {
        acc += parse_float(str[0]);
        acc += parse_float(str[1]);
        acc += parse_float(str[2]);
        acc += parse_float(str[3]);
    }
    UBENCH_DO_NOTHING(&acc);
}

// The C library on the same input, for reference
UBENCH_EX(str, strtod) {
    const char* str[] = {
        "1928123.2767",
        "19.2    ",
        "12323   ",
        "0.000000",
    };

    int64_t num_bytes = 0;
    for (int i = 0; i < (int)ARRAY_SIZE(str); ++i) {
        num_bytes += strlen(str[i]);
    }
    UBENCH_SET_BYTES(num_bytes);

    double acc = 0;
    UBENCH_DO_BENCHMARK() {
        acc += strtod(str[0], NULL);
        acc += strtod(str[1], NULL);
        acc += strtod(str[2], NULL);
        acc += strtod(str[3], NULL);
    }
    UBENCH_DO_NOTHING(&acc);
}
