#include "ubench.h"

#include <core/md_str.h>
#include <core/md_parse.h>
#include <core/md_reader.h>
#include <core/md_file.h>
#include <core/md_allocator.h>
#include <core/md_log.h>

#include <inttypes.h>
#include <stdlib.h>
#include <string.h>

// Lines through a file reader (md_reader.h): the search for the line ends, what md_reader_extract_line costs per byte
UBENCH_EX(str, reader_lines) {
    md_str_t path = MD_STR_INIT(MD_BENCHMARK_DATA_DIR "/centered.gro");
    md_file_reader_t fr;
    if (!md_file_reader_open(&fr, path, MD_MEGABYTES(1), md_get_heap_allocator())) {
        MD_LOG_ERROR("Could not open file '%.*s'", (int)path.len, path.ptr);
        return;
    }

    UBENCH_SET_BYTES(md_reader_size(&fr.r));

    UBENCH_DO_BENCHMARK() {
        md_reader_rewind(&fr.r);
        md_str_t line;
        while (md_reader_extract_line(&fr.r, &line)) {
            // do nothing
        }
    }

    md_file_reader_close(&fr);
}

UBENCH_EX(str, parse_int) {
    md_str_t str[] = {
        MD_STR_INIT("1928123123123"),
        MD_STR_INIT("1123    "),
        MD_STR_INIT("19228123"),
        MD_STR_INIT("1921238123"),
    };

    int64_t num_bytes = 0;
    for (int i = 0; i < (int)MD_COUNTOF(str); ++i) {
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
    for (int i = 0; i < (int)MD_COUNTOF(str); ++i) {
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
    md_str_t str[] = {
        MD_STR_INIT("1928123.2767"),
        MD_STR_INIT("19.2    "),
        MD_STR_INIT("12323   "),
        MD_STR_INIT("0.000000"),
    };

    int64_t num_bytes = 0;
    for (int i = 0; i < (int)MD_COUNTOF(str); ++i) {
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
    for (int i = 0; i < (int)MD_COUNTOF(str); ++i) {
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
