#include "ubench.h"

#include <core/md_str.h>
#include <core/md_parse.h>
#include <core/md_os.h>
#include <core/md_allocator.h>

#include <inttypes.h>

UBENCH_EX(str, read_lines) {
    md_file_o* file = md_file_open(STR(MD_BENCHMARK_DATA_DIR "/centered.gro"), MD_FILE_READ | MD_FILE_BINARY);

    UBENCH_SET_BYTES(md_file_size(file));

    const int64_t cap = MEGABYTES(1);
    char* buf = md_alloc(md_heap_allocator, cap);

    UBENCH_DO_BENCHMARK() {
        md_file_seek(file, 0, MD_FILE_BEG);
        while (md_file_read_lines(file, buf, cap)) {
            // do nothing
        }
    }

    md_free(md_heap_allocator, buf, cap);
}

UBENCH_EX(str, buffered_reader) {
    md_file_o* file = md_file_open(STR(MD_BENCHMARK_DATA_DIR "/centered.gro"), MD_FILE_READ | MD_FILE_BINARY);
    const int64_t cap = MEGABYTES(1);
    char* buf = md_alloc(md_heap_allocator, cap);
    
    md_buffered_reader_t reader = md_buffered_reader_from_file(buf, cap, file);

    UBENCH_SET_BYTES(md_file_size(file));

    UBENCH_DO_BENCHMARK() {
        md_file_seek(file, 0, MD_FILE_BEG);
        str_t line;
        while (md_buffered_reader_extract_line(&line, &reader)) {
            // do nothing
        }
    }

    md_free(md_heap_allocator, buf, cap);
    md_file_close(file);
}

UBENCH_EX(str, parse_int) {
    str_t str[] = {
        STR("1928123123123"),
        STR("1123"),
        STR("19228123"),
        STR("1921238123"),
    };

    int64_t num_bytes = 0;
    for (int i = 0; i < ARRAY_SIZE(str); ++i) {
        num_bytes += str[i].len;
    }
    UBENCH_SET_BYTES(num_bytes);

    uint64_t num_iter = 0;
    uint64_t min = ~0;
    uint64_t max = 0;
    uint64_t avg = 0;

    int64_t acc = 0;
    UBENCH_DO_BENCHMARK() {
        const uint64_t st = __rdtsc();
		acc += parse_int(str[0]);
        acc += parse_int(str[1]);
        acc += parse_int(str[2]);
        acc += parse_int(str[3]);
        const uint64_t et = __rdtsc() - st;
        
        avg += et;
        min = MIN(min, et);
        max = MAX(max, et);
        num_iter += 4;
    }

    printf("avg cycles: %f, min cycles: %"PRIu64", acc: %"PRIi64"\n", (double)(avg) / (double)num_iter, min / 4, acc);
}

UBENCH_EX(str, parse_int_simd) {
    str_t str[] = {
        STR("19312312"),
        STR("1123"),
        STR("19228123"),
        STR("19212381"),
    };

    int64_t num_bytes = 0;
    for (int i = 0; i < ARRAY_SIZE(str); ++i) {
        num_bytes += str[i].len;
    }
    UBENCH_SET_BYTES(num_bytes);

    uint64_t num_iter = 0;
    uint64_t min = ~0;
    uint64_t max = 0;
    uint64_t avg = 0;

    int64_t acc = 0;
    UBENCH_DO_BENCHMARK() {
        const uint64_t st = __rdtsc();
        acc += parse_u32(str[0].ptr, str[0].len);
        acc += parse_u32(str[1].ptr, str[1].len);
        acc += parse_u32(str[2].ptr, str[2].len);
        acc += parse_u32(str[3].ptr, str[3].len);
        const uint64_t et = __rdtsc() - st;

        avg += et;
        min = MIN(min, et);
        max = MAX(max, et);
        num_iter += 4;
    }

    printf("avg cycles: %f, min cycles: %"PRIu64", acc: %"PRIi64"\n", (double)(avg) / (double)num_iter, min / 4, acc);
}

UBENCH_EX(str, parse_float) {
    str_t str[] = {
        STR("1928123.2767"),
        STR("19.2"),
        STR("12323"),
        STR("0.000000"),
    };

    int64_t num_bytes = 0;
    for (int i = 0; i < (int)ARRAY_SIZE(str); ++i) {
        num_bytes += str[i].len;
    }
    UBENCH_SET_BYTES(num_bytes);
    
    uint64_t num_iter = 0;
    uint64_t min = ~0;
    uint64_t max = 0;
    uint64_t avg = 0;

    double acc = 0;
    UBENCH_DO_BENCHMARK() {
        const uint64_t st = __rdtsc();
        acc += parse_float(str[0]);
        acc += parse_float(str[1]);
        acc += parse_float(str[2]);
        acc += parse_float(str[3]);
        const uint64_t et = __rdtsc() - st;

        avg += et;
        min = MIN(min, et);
        max = MAX(max, et);
        num_iter += 4;
    }

    printf("avg cycles: %f, min cycles: %"PRIu64", acc: %f\n", (double)(avg) / (double)num_iter, min / 4, acc);
}

UBENCH_EX(str, parse_float_simd) {
    str_t str[] = {
        STR("1928123.2767"),
        STR("19.2"),
        STR("12323"),
        STR("0.000000"),
    };

    int64_t num_bytes = 0;
    for (int i = 0; i < (int)ARRAY_SIZE(str); ++i) {
        num_bytes += str[i].len;
    }
    UBENCH_SET_BYTES(num_bytes);

    uint64_t num_iter = 0;
    uint64_t min = ~0;
    uint64_t max = 0;
    uint64_t avg = 0;

    double acc = 0;
    UBENCH_DO_BENCHMARK() {
        const uint64_t st = __rdtsc();
        acc += parse_float_wide(str[0].ptr, str[0].len);
        acc += parse_float_wide(str[1].ptr, str[1].len);
        acc += parse_float_wide(str[2].ptr, str[2].len);
        acc += parse_float_wide(str[3].ptr, str[3].len);
        const uint64_t et = __rdtsc() - st;

        avg += et;
        min = MIN(min, et);
        max = MAX(max, et);
        num_iter += 4;
    }

    printf("avg cycles: %f, min cycles: %"PRIu64", acc: %f\n", (double)(avg) / (double)num_iter, min / 4, acc);
}