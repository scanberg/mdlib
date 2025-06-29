#include "ubench.h"

#include <core/md_str.h>
#include <core/md_parse.h>
#include <core/md_os.h>
#include <core/md_allocator.h>
#include <core/md_log.h>

#include <inttypes.h>

UBENCH_EX(str, buffered_reader) {
    str_t path = STR_LIT(MD_BENCHMARK_DATA_DIR "/centered.gro");
    md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);
    if (!file) {
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
    md_file_close(file);
}

UBENCH_EX(str, parse_int) {
    str_t str[] = {
        STR_LIT("1928123123123"),
        STR_LIT("1123    "),
        STR_LIT("19228123"),
        STR_LIT("1921238123"),
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

UBENCH_EX(str, parse_int_simd) {
    str_t str[] = {
        STR_LIT("19312312"),
        STR_LIT("1123    "),
        STR_LIT("19228123"),
        STR_LIT("19212381"),
    };

    int64_t num_bytes = 0;
    for (int i = 0; i < (int)ARRAY_SIZE(str); ++i) {
        num_bytes += str[i].len;
    }
    UBENCH_SET_BYTES(num_bytes);

    size_t acc = 0;
    UBENCH_DO_BENCHMARK() {
        acc += parse_u32(str[0].ptr, str[0].len);
        acc += parse_u32(str[1].ptr, str[1].len);
        acc += parse_u32(str[2].ptr, str[2].len);
        acc += parse_u32(str[3].ptr, str[3].len);
    }
    UBENCH_DO_NOTHING(&acc);
}

UBENCH_EX(str, parse_float) {
    str_t str[] = {
        STR_LIT("1928123.2767"),
        STR_LIT("19.2    "),
        STR_LIT("12323   "),
        STR_LIT("0.000000"),
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

UBENCH_EX(str, parse_float_simd) {
    str_t str[] = {
        STR_LIT("1928123.2767    "),
        STR_LIT("19.2            "),
        STR_LIT("12323           "),
        STR_LIT("0.000000        "),
    };

    int64_t num_bytes = 0;
    for (int i = 0; i < (int)ARRAY_SIZE(str); ++i) {
        num_bytes += str[i].len;
    }
    UBENCH_SET_BYTES(num_bytes);

    double acc = 0;
    UBENCH_DO_BENCHMARK() {
        acc += parse_float_wide(str[0].ptr, str[0].len);
        acc += parse_float_wide(str[1].ptr, str[1].len);
        acc += parse_float_wide(str[2].ptr, str[2].len);
        acc += parse_float_wide(str[3].ptr, str[3].len);
    }
    UBENCH_DO_NOTHING(&acc);
}
